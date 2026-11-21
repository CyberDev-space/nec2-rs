//! Structure geometry: wire segments, surface patches and the way they
//! connect.
//!
//! [`Geometry`] collects segments and patches the way the geometry section of
//! a NEC deck does (GW, GA, GH, SP, SM, SC, GM, GX, GR, GS).
//! [`Geometry::connect`] ends the input (GE): it finds the segment ends that
//! touch each other or a patch, makes them coincide exactly and returns an
//! immutable [`Structure`].
//!
//! Ported from `geometry.c` of nec2c 1.3.1 (`wire`, `arc`, `helix`, `patch`,
//! `subph`, `move`, `reflc`, `conect`, `isegno` and the GE part of `datagn`).

use crate::consts::{PCHCON, PI, SMIN, TA, TP};
use crate::error::{Error, Result};

/// A surface patch: a small flat piece of a closed conducting surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Patch {
    /// Centre.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// Two unit vectors in the surface, at right angles to each other.
    pub t1x: f64,
    pub t1y: f64,
    pub t1z: f64,
    pub t2x: f64,
    pub t2y: f64,
    pub t2z: f64,
    /// Area.
    pub area: f64,
    /// +1 or -1: the outward normal is `salp` times t1 x t2.
    pub salp: f64,
}

impl Patch {
    /// Unit normal that points out of the surface.
    pub fn normal(&self) -> [f64; 3] {
        [
            (self.t1y * self.t2z - self.t1z * self.t2y) * self.salp,
            (self.t1z * self.t2x - self.t1x * self.t2z) * self.salp,
            (self.t1x * self.t2y - self.t1y * self.t2x) * self.salp,
        ]
    }
}

/// Divides patch number `nx` (from 1) into four patches around its centre,
/// where a wire is connected (`subph`). The four take its place in the list.
fn subph(patches: &mut Vec<Patch>, mp: &mut usize, nx: usize) {
    let mi = nx - 1;
    let whole = patches[mi];
    let xa = whole.area / 4.0;
    let xst = xa.sqrt() / 2.0;
    let (mut xt, mut yt) = (xst, xst);
    let mut quarters = [whole; 4];
    for (ix, quarter) in quarters.iter_mut().enumerate() {
        quarter.x = whole.x + xt * whole.t1x + yt * whole.t2x;
        quarter.y = whole.y + xt * whole.t1y + yt * whole.t2y;
        quarter.z = whole.z + xt * whole.t1z + yt * whole.t2z;
        quarter.area = xa;
        if ix == 1 {
            yt = -yt;
        }
        if ix == 0 || ix == 2 {
            xt = -xt;
        }
    }
    patches.splice(mi..=mi, quarters);
    if nx <= *mp {
        *mp += 3;
    }
}

/// Wire segments and surface patches under construction.
///
/// Dimensions are metres (or any unit, if [`Geometry::scale`] converts them to
/// metres before [`Geometry::connect`]).
#[derive(Debug, Clone, Default)]
pub struct Geometry {
    /// Number of segments and of patches in one symmetric cell.
    np: usize,
    mp: usize,
    /// Symmetry flag: 0 none, 1 planes without Z, 2 planes with Z, -1 rotation.
    ipsym: i32,
    patches: Vec<Patch>,
    itag: Vec<i32>,
    x1: Vec<f64>,
    y1: Vec<f64>,
    z1: Vec<f64>,
    x2: Vec<f64>,
    y2: Vec<f64>,
    z2: Vec<f64>,
    bi: Vec<f64>,
}

/// What to do when the centre of a segment lies inside the volume of another
/// wire. NEC-2 gives meaningless results for such structures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OverlapCheck {
    /// Refuse the structure.
    Error,
    /// Accept the structure, as NEC-2 does, and list the pairs in
    /// [`Structure::warnings`].
    #[default]
    Warn,
    /// Do not test.
    Off,
}

/// A segment whose centre lies inside the volume of another wire.
#[derive(Debug, Clone, PartialEq)]
pub struct Overlap {
    /// The segment whose centre is inside (1-based) and its tag.
    pub segment: usize,
    pub tag: i32,
    /// The segment whose volume it is in (1-based) and its tag.
    pub inside: usize,
    pub inside_tag: i32,
    /// Distance of the centre from the axis of the other segment.
    pub distance: f64,
    /// Radius of the other segment.
    pub radius: f64,
}

impl std::fmt::Display for Overlap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SEGMENT {} (TAG {}) HAS ITS CENTRE INSIDE THE VOLUME OF SEGMENT {} (TAG {}): DISTANCE {:.4E} IS LESS THAN THE RADIUS {:.4E}",
            self.segment, self.tag, self.inside, self.inside_tag, self.distance, self.radius
        )
    }
}

/// End coordinates of all segments.
struct Ends<'a> {
    x1: &'a [f64],
    y1: &'a [f64],
    z1: &'a [f64],
    x2: &'a [f64],
    y2: &'a [f64],
    z2: &'a [f64],
}

impl Ends<'_> {
    /// End `k` of the structure: `2 i` is end 1 of segment `i`, `2 i + 1` its end 2.
    fn point(&self, k: usize) -> [f64; 3] {
        let i = k / 2;
        if k % 2 == 0 { [self.x1[i], self.y1[i], self.z1[i]] } else { [self.x2[i], self.y2[i], self.z2[i]] }
    }
}

/// Groups the segment ends into nodes. Two ends belong to one node if the
/// sum of their coordinate differences is at most the larger of the
/// connection tolerances of their segments; the grouping is transitive.
/// Returns the node of every end.
fn end_nodes(ends: &Ends, si: &[f64]) -> Vec<usize> {
    let count = 2 * si.len();
    let tolerance = |k: usize| si[k / 2] * SMIN;
    let largest = si.iter().fold(0.0f64, |a, &b| a.max(b)) * SMIN;
    let mut parent: Vec<usize> = (0..count).collect();
    fn find(parent: &mut [usize], mut k: usize) -> usize {
        while parent[k] != k {
            parent[k] = parent[parent[k]];
            k = parent[k];
        }
        k
    }

    // Only ends that are close along X can belong together.
    let mut order: Vec<usize> = (0..count).collect();
    order.sort_by(|&a, &b| ends.point(a)[0].total_cmp(&ends.point(b)[0]));
    for (position, &a) in order.iter().enumerate() {
        let pa = ends.point(a);
        for &b in &order[position + 1..] {
            let pb = ends.point(b);
            if pb[0] - pa[0] > largest {
                break;
            }
            let separation = (pa[0] - pb[0]).abs() + (pa[1] - pb[1]).abs() + (pa[2] - pb[2]).abs();
            if separation <= tolerance(a).max(tolerance(b)) {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                parent[ra] = rb;
            }
        }
    }
    (0..count).map(|k| find(&mut parent, k)).collect()
}

/// Distance of `point` from the segment between `a` and `b`.
fn distance_to_axis(point: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let p = [point[0] - a[0], point[1] - a[1], point[2] - a[2]];
    let length2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    let t = ((p[0] * d[0] + p[1] * d[1] + p[2] * d[2]) / length2).clamp(0.0, 1.0);
    let r = [p[0] - t * d[0], p[1] - t * d[1], p[2] - t * d[2]];
    (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt()
}

/// Finds the segments whose centre lies inside the volume of another
/// segment. Two segments that meet at exactly one node form a junction and
/// are not tested against each other.
fn overlaps(ends: &Ends, centre: [&[f64]; 3], si: &[f64], bi: &[f64], itag: &[i32]) -> Vec<Overlap> {
    let n = si.len();
    let nodes = end_nodes(ends, si);
    let centre_of = |i: usize| [centre[0][i], centre[1][i], centre[2][i]];
    let mut found = Vec::new();
    for a in 0..n {
        let ca = centre_of(a);
        for b in a + 1..n {
            let cb = centre_of(b);
            // Segments further apart than their reach cannot overlap.
            let reach = (0.5 * si[a] + bi[a]).max(0.5 * si[b] + bi[b]);
            let apart = [(ca[0] - cb[0]).abs(), (ca[1] - cb[1]).abs(), (ca[2] - cb[2]).abs()];
            if apart.iter().any(|&d| d > reach) {
                continue;
            }
            let (a1, a2, b1, b2) = (nodes[2 * a], nodes[2 * a + 1], nodes[2 * b], nodes[2 * b + 1]);
            // The number of nodes the two segments have in common.
            let first = usize::from(a1 == b1 || a1 == b2);
            let shared = if a1 == a2 { first } else { first + usize::from(a2 == b1 || a2 == b2) };
            if shared == 1 {
                continue;
            }
            // One entry for the pair, whichever centre is inside.
            let hit = [(a, b, ca), (b, a, cb)].into_iter().find_map(|(inner, outer, point)| {
                let distance = distance_to_axis(point, ends.point(2 * outer), ends.point(2 * outer + 1));
                (distance < bi[outer]).then(|| Overlap {
                    segment: inner + 1,
                    tag: itag[inner],
                    inside: outer + 1,
                    inside_tag: itag[outer],
                    distance,
                    radius: bi[outer],
                })
            });
            found.extend(hit);
        }
    }
    found
}

/// The pairs as lines of text. A long list is cut short.
pub(crate) fn overlap_lines(pairs: &[Overlap]) -> Vec<String> {
    const SHOWN: usize = 10;
    let mut lines: Vec<String> = pairs.iter().take(SHOWN).map(Overlap::to_string).collect();
    if pairs.len() > SHOWN {
        lines.push(format!("AND {} MORE PAIRS", pairs.len() - SHOWN));
    }
    lines
}

/// What [`Geometry::helix`] reports about the wire it generated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HelixInfo {
    /// The radius changes along the axis; `cone_angle` is in radians.
    Spiral { cone_angle: f64 },
    /// Constant radius; `pitch_angle` is in degrees, `turn_length` in input units.
    Helix { pitch_angle: f64, turn_length: f64 },
}

/// Returns the 1-based number of the `mx`-th segment that carries tag `itagi`.
/// With `itagi == 0` the result is segment number `mx`.
pub fn isegno(itag: &[i32], itagi: i32, mx: i32) -> Result<usize> {
    if mx <= 0 {
        return Err(Error::ZeroSegmentPosition);
    }
    if itagi == 0 {
        // nec2c returns `mx` unchecked here; a number past the last segment
        // would index outside the structure.
        return if mx as usize <= itag.len() { Ok(mx as usize) } else { Err(Error::SegmentOutOfRange(mx)) };
    }
    let mut icnt = 0;
    for (i, &tag) in itag.iter().enumerate() {
        if tag == itagi {
            icnt += 1;
            if icnt == mx {
                return Ok(i + 1);
            }
        }
    }
    Err(Error::NoSegmentWithTag(itagi))
}

impl Geometry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of segments entered so far.
    pub fn segment_count(&self) -> usize {
        self.bi.len()
    }

    /// Number of patches entered so far.
    pub fn patch_count(&self) -> usize {
        self.patches.len()
    }

    /// New patches (SP, SM and SC), from `patch` of nec2c.
    ///
    /// With `nx == 0`, `ny` selects the shape of one patch. 1: any shape,
    /// with the centre `c1`, the elevation and the azimuth of the outward
    /// normal in radians and the area in `c2`. 2: a rectangle with the corners
    /// `c1`, `c2`, `c3`. 3: a triangle with these corners. 4: a quadrilateral
    /// with the corners `c1` to `c4`. With `nx > 0` a rectangle with the
    /// corners `c1`, `c2`, `c3` is covered by `nx` times `ny` patches.
    pub fn patch(&mut self, nx: i32, ny: i32, c1: [f64; 3], c2: [f64; 3], c3: [f64; 3], c4: [f64; 3]) -> Result<()> {
        let [ax1, ay1, az1] = c1;
        let [ax2, ay2, az2] = c2;
        let [ax3, ay3, az3] = c3;
        let [ax4, ay4, az4] = c4;
        let ntp = if nx > 0 { 2 } else { ny };
        let (mut s1x, mut s1y, mut s1z, mut s2x, mut s2y, mut s2z) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let (xnv, ynv, znv, xa);
        let (px, py, pz, pbi, t1x, t1y, t1z);

        if ntp <= 1 {
            px = ax1;
            py = ay1;
            pz = az1;
            pbi = az2;
            let cos_elevation = ax2.cos();
            xnv = cos_elevation * ay2.cos();
            ynv = cos_elevation * ay2.sin();
            znv = ax2.sin();
            xa = (xnv * xnv + ynv * ynv).sqrt();
            if xa >= 1.0e-6 {
                t1x = -ynv / xa;
                t1y = xnv / xa;
                t1z = 0.0;
            } else {
                t1x = 1.0;
                t1y = 0.0;
                t1z = 0.0;
            }
        } else {
            s1x = ax2 - ax1;
            s1y = ay2 - ay1;
            s1z = az2 - az1;
            s2x = ax3 - ax2;
            s2y = ay3 - ay2;
            s2z = az3 - az2;
            if nx != 0 {
                s1x /= f64::from(nx);
                s1y /= f64::from(nx);
                s1z /= f64::from(nx);
                s2x /= f64::from(ny);
                s2y /= f64::from(ny);
                s2z /= f64::from(ny);
            }
            let (xn, yn, zn) = (s1y * s2z - s1z * s2y, s1z * s2x - s1x * s2z, s1x * s2y - s1y * s2x);
            xa = (xn * xn + yn * yn + zn * zn).sqrt();
            xnv = xn / xa;
            ynv = yn / xa;
            znv = zn / xa;
            let xst = (s1x * s1x + s1y * s1y + s1z * s1z).sqrt();
            t1x = s1x / xst;
            t1y = s1y / xst;
            t1z = s1z / xst;

            if ntp <= 2 {
                px = ax1 + 0.5 * (s1x + s2x);
                py = ay1 + 0.5 * (s1y + s2y);
                pz = az1 + 0.5 * (s1z + s2z);
                pbi = xa;
            } else if ntp != 4 {
                px = (ax1 + ax2 + ax3) / 3.0;
                py = (ay1 + ay2 + ay3) / 3.0;
                pz = (az1 + az2 + az3) / 3.0;
                pbi = 0.5 * xa;
            } else {
                // Two triangles; their normals must agree.
                let (d1x, d1y, d1z) = (ax3 - ax1, ay3 - ay1, az3 - az1);
                let (d2x, d2y, d2z) = (ax4 - ax1, ay4 - ay1, az4 - az1);
                let xn2 = d1y * d2z - d1z * d2y;
                let yn2 = d1z * d2x - d1x * d2z;
                let zn2 = d1x * d2y - d1y * d2x;
                let xst = (xn2 * xn2 + yn2 * yn2 + zn2 * zn2).sqrt();
                let salpn = 1.0 / (3.0 * (xa + xst));
                px = (xa * (ax1 + ax2 + ax3) + xst * (ax1 + ax3 + ax4)) * salpn;
                py = (xa * (ay1 + ay2 + ay3) + xst * (ay1 + ay3 + ay4)) * salpn;
                pz = (xa * (az1 + az2 + az3) + xst * (az1 + az3 + az4)) * salpn;
                pbi = 0.5 * (xa + xst);
                let cosine = (xnv * xn2 + ynv * yn2 + znv * zn2) / xst;
                if cosine <= 0.9998 {
                    return Err(Error::PatchNotPlanar);
                }
            }
        }

        let first = Patch {
            x: px,
            y: py,
            z: pz,
            t1x,
            t1y,
            t1z,
            t2x: ynv * t1z - znv * t1y,
            t2y: znv * t1x - xnv * t1z,
            t2z: xnv * t1y - ynv * t1x,
            area: pbi,
            salp: 1.0,
        };

        if nx == 0 {
            self.patches.push(first);
        } else {
            let mut xn2 = first.x - s1x - s2x;
            let mut yn2 = first.y - s1y - s2y;
            let mut zn2 = first.z - s1z - s2z;
            for _ in 0..ny {
                xn2 += s2x;
                yn2 += s2y;
                zn2 += s2z;
                for ix in 1..=nx {
                    let xst = f64::from(ix);
                    self.patches.push(Patch {
                        x: xn2 + xst * s1x,
                        y: yn2 + xst * s1y,
                        z: zn2 + xst * s1z,
                        area: xa,
                        ..first
                    });
                }
            }
        }
        self.reset_symmetry();
        Ok(())
    }

    fn push(&mut self, tag: i32, end1: [f64; 3], end2: [f64; 3], rad: f64) {
        self.itag.push(tag);
        self.x1.push(end1[0]);
        self.y1.push(end1[1]);
        self.z1.push(end1[2]);
        self.x2.push(end2[0]);
        self.y2.push(end2[1]);
        self.z2.push(end2[2]);
        self.bi.push(rad);
    }

    /// New segments end any symmetry declared so far.
    fn reset_symmetry(&mut self) {
        self.np = self.segment_count();
        self.mp = self.patches.len();
        self.ipsym = 0;
    }

    /// Straight wire of `ns` segments (GW, with the taper of GC).
    ///
    /// `rdel` is the ratio of the lengths of adjacent segments and `rrad` the
    /// ratio of their radii; both are 1 for a uniform wire. A wire of zero
    /// segments is ignored, as in NEC-2.
    #[allow(clippy::too_many_arguments)]
    pub fn wire(
        &mut self,
        tag: i32,
        ns: i32,
        end1: [f64; 3],
        end2: [f64; 3],
        rad: f64,
        rdel: f64,
        rrad: f64,
    ) -> Result<()> {
        if ns < 0 {
            return Err(Error::InvalidSegmentCount(ns));
        }
        if ns == 0 {
            self.reset_symmetry();
            return Ok(());
        }
        let [xw1, yw1, zw1] = end1;
        let [xw2, yw2, zw2] = end2;
        let mut xd = xw2 - xw1;
        let mut yd = yw2 - yw1;
        let mut zd = zw2 - zw1;
        let mut delz;
        let rd;
        if (rdel - 1.0).abs() >= 1.0e-6 {
            delz = (xd * xd + yd * yd + zd * zd).sqrt();
            xd /= delz;
            yd /= delz;
            zd /= delz;
            delz = delz * (1.0 - rdel) / (1.0 - rdel.powf(f64::from(ns)));
            rd = rdel;
        } else {
            let fns = f64::from(ns);
            xd /= fns;
            yd /= fns;
            zd /= fns;
            delz = 1.0;
            rd = 1.0;
        }

        let mut radz = rad;
        let (mut xs1, mut ys1, mut zs1) = (xw1, yw1, zw1);
        for _ in 0..ns {
            let xs2 = xs1 + xd * delz;
            let ys2 = ys1 + yd * delz;
            let zs2 = zs1 + zd * delz;
            self.push(tag, [xs1, ys1, zs1], [xs2, ys2, zs2], radz);
            delz *= rd;
            radz *= rrad;
            xs1 = xs2;
            ys1 = ys2;
            zs1 = zs2;
        }
        if ns >= 1 {
            let last = self.segment_count() - 1;
            self.x2[last] = xw2;
            self.y2[last] = yw2;
            self.z2[last] = zw2;
        }
        self.reset_symmetry();
        Ok(())
    }

    /// Circular arc of `ns` segments in the XZ plane, centred on the origin (GA).
    /// Angles are degrees from the X axis.
    pub fn arc(&mut self, tag: i32, ns: i32, rada: f64, ang1: f64, ang2: f64, rad: f64) -> Result<()> {
        if ns < 0 {
            return Err(Error::InvalidSegmentCount(ns));
        }
        if ns >= 1 {
            if (ang2 - ang1).abs() >= 360.00001 {
                return Err(Error::ArcAngleExceeds360);
            }
            let mut ang = ang1 * TA;
            let dang = (ang2 - ang1) * TA / f64::from(ns);
            let mut xs1 = rada * ang.cos();
            let mut zs1 = rada * ang.sin();
            for _ in 0..ns {
                ang += dang;
                let xs2 = rada * ang.cos();
                let zs2 = rada * ang.sin();
                self.push(tag, [xs1, 0.0, zs1], [xs2, 0.0, zs2], rad);
                xs1 = xs2;
                zs1 = zs2;
            }
        }
        self.reset_symmetry();
        Ok(())
    }

    /// Helix or spiral of `ns` segments along the Z axis (GH).
    ///
    /// `s` is the spacing between turns, `hl` the total length (negative for a
    /// left-handed helix), `a1`/`b1` the radii along X and Y at z = 0 and
    /// `a2`/`b2` the radii at z = `hl`.
    #[allow(clippy::too_many_arguments)]
    pub fn helix(
        &mut self,
        tag: i32,
        ns: i32,
        s: f64,
        hl: f64,
        a1: f64,
        b1: f64,
        a2: f64,
        b2: f64,
        rad: f64,
    ) -> Result<Option<HelixInfo>> {
        if ns < 0 {
            return Err(Error::InvalidSegmentCount(ns));
        }
        if ns == 0 {
            self.reset_symmetry();
            return Ok(None);
        }
        if s == 0.0 {
            return Err(Error::ZeroHelixSpacing);
        }
        if hl == 0.0 {
            // Flat spiral in the plane z = 0: `s` is the number of segments
            // per turn and the radii change linearly from node to node.
            let dphi = 2.0 * std::f64::consts::PI / s;
            let node = |k: i32| {
                let t = f64::from(k) / f64::from(ns);
                let phi = f64::from(k) * dphi;
                [(a1 + (a2 - a1) * t) * phi.cos(), (b1 + (b2 - b1) * t) * phi.sin(), 0.0]
            };
            for k in 1..=ns {
                self.push(tag, node(k - 1), node(k), rad);
            }
            self.reset_symmetry();
            return Ok(None);
        }
        let (mut b1, mut b2) = (b1, b2);
        let zinc = (hl / f64::from(ns)).abs();
        let mut z1 = 0.0;
        for i in 0..ns {
            if i != 0 {
                z1 += zinc;
            }
            let z2 = z1 + zinc;
            let (mut x1, mut y1, mut x2, mut y2);
            if a2 == a1 {
                if b1 == 0.0 {
                    b1 = a1;
                }
                x1 = a1 * (2.0 * PI * z1 / s).cos();
                y1 = b1 * (2.0 * PI * z1 / s).sin();
                x2 = a1 * (2.0 * PI * z2 / s).cos();
                y2 = b1 * (2.0 * PI * z2 / s).sin();
            } else {
                if b2 == 0.0 {
                    b2 = a2;
                }
                x1 = (a1 + (a2 - a1) * z1 / hl.abs()) * (2.0 * PI * z1 / s).cos();
                y1 = (b1 + (b2 - b1) * z1 / hl.abs()) * (2.0 * PI * z1 / s).sin();
                x2 = (a1 + (a2 - a1) * z2 / hl.abs()) * (2.0 * PI * z2 / s).cos();
                y2 = (b1 + (b2 - b1) * z2 / hl.abs()) * (2.0 * PI * z2 / s).sin();
            }
            if hl <= 0.0 {
                std::mem::swap(&mut x1, &mut y1);
                std::mem::swap(&mut x2, &mut y2);
            }
            self.push(tag, [x1, y1, z1], [x2, y2, z2], rad);
        }
        self.reset_symmetry();

        if a2 != a1 {
            // The angle whose tangent is the change of the radius over the
            // length. The formula of NEC-2 says the same, but for a spiral
            // that closes to a point it divides zero by zero.
            let cone_angle = if a2 == 0.0 {
                (-a1 / hl.abs()).atan()
            } else {
                (a2 / (hl.abs() + (hl.abs() * a1) / (a2 - a1))).atan()
            };
            return Ok(Some(HelixInfo::Spiral { cone_angle }));
        }
        let (pitch_angle, turn_length);
        if a1 == b1 {
            let hdia = 2.0 * a1;
            let pitch = (s / (PI * hdia)).atan();
            turn_length = hdia * PI / pitch.cos();
            pitch_angle = 180.0 * pitch / PI;
        } else {
            let (hmaj, hmin) = if a1 >= b1 { (2.0 * a1, 2.0 * b1) } else { (2.0 * b1, 2.0 * a1) };
            // Printed for information only; the expression is the one nec2c uses.
            let hdia = ((hmaj * hmaj + hmin * hmin) / 2.0 * hmaj).sqrt();
            turn_length = 2.0 * PI * hdia;
            pitch_angle = (180.0 / PI) * (s / (PI * hdia)).atan();
        }
        Ok(Some(HelixInfo::Helix { pitch_angle, turn_length }))
    }

    /// Rotates the structure about the X, Y and Z axes by `rox`, `roy`, `roz`
    /// radians and then shifts it by `(xs, ys, zs)` (GM).
    ///
    /// With `nrpt == 0` the segments from the first one tagged `its` onward
    /// (all segments if `its == 0`) are moved in place. Otherwise `nrpt` new
    /// copies are appended, each one moved relative to the previous copy.
    /// Non-zero tags of moved or copied segments are increased by `itgi`.
    #[allow(clippy::too_many_arguments)]
    pub fn move_structure(
        &mut self,
        rox: f64,
        roy: f64,
        roz: f64,
        xs: f64,
        ys: f64,
        zs: f64,
        its: i32,
        nrpt: i32,
        itgi: i32,
    ) -> Result<()> {
        if nrpt < 0 {
            return Err(Error::NegativeField("REPETITION COUNT", nrpt));
        }
        if rox.abs() + roy.abs() > 1.0e-10 {
            self.ipsym *= 3;
        }
        let (sps, cps) = (rox.sin(), rox.cos());
        let (sth, cth) = (roy.sin(), roy.cos());
        let (sph, cph) = (roz.sin(), roz.cos());
        let xx = cph * cth;
        let xy = cph * sth * sps - sph * cps;
        let xz = cph * sth * cps + sph * sps;
        let yx = sph * cth;
        let yy = sph * sth * sps + cph * cps;
        let yz = sph * sth * cps - cph * sps;
        let zx = -sth;
        let zy = cth * sps;
        let zz = cth * cps;
        let moved = |p: [f64; 3]| {
            let [xi, yi, zi] = p;
            [xi * xx + yi * xy + zi * xz + xs, xi * yx + yi * yy + zi * yz + ys, xi * zx + yi * zy + zi * zz + zs]
        };

        let mut ix = 1;
        if self.segment_count() > 0 {
            let mut i1 = isegno(&self.itag, its, 1)?;
            ix = i1;
            let repeats = if nrpt == 0 { 1 } else { nrpt.max(0) };
            for _ in 0..repeats {
                let n = self.segment_count();
                for i in (i1 - 1)..n {
                    let end1 = moved([self.x1[i], self.y1[i], self.z1[i]]);
                    let end2 = moved([self.x2[i], self.y2[i], self.z2[i]]);
                    let tag = if self.itag[i] != 0 { self.itag[i] + itgi } else { 0 };
                    if nrpt == 0 {
                        self.itag[i] = tag;
                        [self.x1[i], self.y1[i], self.z1[i]] = end1;
                        [self.x2[i], self.y2[i], self.z2[i]] = end2;
                    } else {
                        self.push(tag, end1, end2, self.bi[i]);
                    }
                }
                i1 = n + 1;
            }
        }

        // All patches move, whatever the first tag.
        if !self.patches.is_empty() {
            let turned = |v: [f64; 3]| {
                let [xi, yi, zi] = v;
                [xi * xx + yi * xy + zi * xz, xi * yx + yi * yy + zi * yz, xi * zx + yi * zy + zi * zz]
            };
            let repeats = if nrpt == 0 { 1 } else { nrpt };
            let mut i1 = 0;
            for _ in 0..repeats {
                let m = self.patches.len();
                for i in i1..m {
                    let old = self.patches[i];
                    let [x, y, z] = moved([old.x, old.y, old.z]);
                    let [t1x, t1y, t1z] = turned([old.t1x, old.t1y, old.t1z]);
                    let [t2x, t2y, t2z] = turned([old.t2x, old.t2y, old.t2z]);
                    let new = Patch { x, y, z, t1x, t1y, t1z, t2x, t2y, t2z, ..old };
                    if nrpt == 0 {
                        self.patches[i] = new;
                    } else {
                        self.patches.push(new);
                    }
                }
                i1 = m;
            }
        }

        if nrpt == 0 && ix == 1 {
            return Ok(());
        }
        self.reset_symmetry();
        Ok(())
    }

    /// Appends the mirror image of the whole structure in the plane normal to
    /// coordinate `axis` (0 = X, 1 = Y, 2 = Z).
    fn mirror(&mut self, axis: usize, iti: i32) -> Result<()> {
        let n = self.segment_count();
        for i in 0..n {
            let mut end1 = [self.x1[i], self.y1[i], self.z1[i]];
            let mut end2 = [self.x2[i], self.y2[i], self.z2[i]];
            let (e1, e2) = (end1[axis], end2[axis]);
            if e1.abs() + e2.abs() <= 1.0e-5 || e1 * e2 < -1.0e-6 {
                return Err(Error::SegmentInSymmetryPlane(i + 1));
            }
            end1[axis] = -e1;
            end2[axis] = -e2;
            let tag = if self.itag[i] != 0 { self.itag[i] + iti } else { 0 };
            self.push(tag, end1, end2, self.bi[i]);
        }

        // The image of a patch: its normal is turned over with it.
        for i in 0..self.patches.len() {
            let mut image = self.patches[i];
            let (centre, t1, t2) = match axis {
                0 => (&mut image.x, &mut image.t1x, &mut image.t2x),
                1 => (&mut image.y, &mut image.t1y, &mut image.t2y),
                _ => (&mut image.z, &mut image.t1z, &mut image.t2z),
            };
            if centre.abs() <= 1.0e-10 {
                return Err(Error::PatchInSymmetryPlane(i + 1));
            }
            *centre = -*centre;
            *t1 = -*t1;
            *t2 = -*t2;
            image.salp = -image.salp;
            self.patches.push(image);
        }
        Ok(())
    }

    /// Completes a structure that has plane symmetry by reflecting it in the
    /// coordinate planes normal to the chosen axes (GX). Tags of the images
    /// are increased by `itx`, doubled after each reflection.
    pub fn reflect(&mut self, along_x: bool, along_y: bool, along_z: bool, itx: i32) -> Result<()> {
        self.reset_symmetry();
        if !(along_x || along_y || along_z) {
            return Ok(());
        }
        // The tag increment doubles only while there are segments.
        let doubling = if self.segment_count() > 0 { 2 } else { 1 };
        let mut iti = itx;
        self.ipsym = 1;
        if along_z {
            self.ipsym = 2;
            self.mirror(2, iti)?;
            iti *= doubling;
        }
        if along_y {
            self.mirror(1, iti)?;
            iti *= doubling;
        }
        if along_x {
            self.mirror(0, iti)?;
        }
        Ok(())
    }

    /// Completes a structure that has rotational symmetry by copying it
    /// `nop - 1` times around the Z axis (GR). Tags grow by `itx` per copy.
    pub fn rotate(&mut self, itx: i32, nop: i32) -> Result<()> {
        if nop < 1 {
            return Err(Error::InvalidRotationCount(nop));
        }
        if nop == 1 {
            // One sector is the structure itself.
            return Ok(());
        }
        self.reset_symmetry();
        self.ipsym = -1;
        let sam = TP / f64::from(nop);
        let (cs, ss) = (sam.cos(), sam.sin());
        let np = self.np;
        for i in np..np * nop as usize {
            let k = i - np;
            let end1 = [self.x1[k] * cs - self.y1[k] * ss, self.x1[k] * ss + self.y1[k] * cs, self.z1[k]];
            let end2 = [self.x2[k] * cs - self.y2[k] * ss, self.x2[k] * ss + self.y2[k] * cs, self.z2[k]];
            let tag = if self.itag[k] != 0 { self.itag[k] + itx } else { 0 };
            self.push(tag, end1, end2, self.bi[k]);
        }
        let mp = self.mp;
        for i in mp..mp * nop as usize {
            let k = self.patches[i - mp];
            self.patches.push(Patch {
                x: k.x * cs - k.y * ss,
                y: k.x * ss + k.y * cs,
                t1x: k.t1x * cs - k.t1y * ss,
                t1y: k.t1x * ss + k.t1y * cs,
                t2x: k.t2x * cs - k.t2y * ss,
                t2y: k.t2x * ss + k.t2y * cs,
                ..k
            });
        }
        Ok(())
    }

    /// Multiplies all dimensions, including wire radii, by `factor` (GS).
    pub fn scale(&mut self, factor: f64) {
        let Self { x1, y1, z1, x2, y2, z2, bi, patches, .. } = self;
        for values in [x1, y1, z1, x2, y2, z2, bi] {
            for v in values.iter_mut() {
                *v *= factor;
            }
        }
        for patch in patches.iter_mut() {
            patch.x *= factor;
            patch.y *= factor;
            patch.z *= factor;
            patch.area *= factor * factor;
        }
    }

    /// Ends the geometry input (GE).
    ///
    /// `ignd` is the ground-plane flag of the GE card: 0 means no ground,
    /// a positive value means that wires ending on z = 0 connect to their
    /// image, a negative value means a ground is present but the current
    /// goes to zero where a wire touches it.
    pub fn connect(self, ignd: i32) -> Result<Structure> {
        self.connect_with(ignd, OverlapCheck::Error)
    }

    /// Like [`Geometry::connect`], with a choice of what happens when wires
    /// occupy the same space.
    pub fn connect_with(self, ignd: i32, overlap: OverlapCheck) -> Result<Structure> {
        let n = self.segment_count();
        if n + self.patches.len() == 0 {
            return Err(Error::EmptyStructure);
        }
        let Geometry {
            mut np,
            mut mp,
            mut ipsym,
            mut patches,
            itag,
            mut x1,
            mut y1,
            mut z1,
            mut x2,
            mut y2,
            mut z2,
            bi,
        } = self;

        // Every segment needs finite ends, a length and a radius. The tests
        // are written so that a NaN fails them.
        for i in 0..n {
            let (dx, dy, dz) = (x2[i] - x1[i], y2[i] - y1[i], z2[i] - z1[i]);
            let length = (dx * dx + dy * dy + dz * dz).sqrt();
            let ends_finite = [x1[i], y1[i], z1[i], x2[i], y2[i], z2[i]].iter().all(|c| c.is_finite());
            let valid = ends_finite && length.is_finite() && length > 1.0e-20 && bi[i].is_finite() && bi[i] > 0.0;
            if !valid {
                return Err(Error::SegmentData { segment: i + 1, tag: itag[i] });
            }
        }
        // Every patch needs a finite centre, finite tangents and a finite
        // area. NEC-2 solves a structure with a patch of no area.
        let mut warnings = Vec::new();
        for (i, p) in patches.iter().enumerate() {
            let values = [p.x, p.y, p.z, p.t1x, p.t1y, p.t1z, p.t2x, p.t2y, p.t2z, p.area];
            if !values.iter().all(|c| c.is_finite()) {
                return Err(Error::PatchData(i + 1));
            }
            if p.area <= 0.0 {
                warnings.push(format!("PATCH {} HAS NO AREA", i + 1));
            }
        }

        if ignd != 0 {
            if ipsym == 2 {
                np *= 2;
                mp *= 2;
            }
            if ipsym.abs() > 2 {
                np = n;
                mp = patches.len();
            }
            if np > n {
                return Err(Error::Symmetry { cell: np, total: n });
            }
            if np == n && mp == patches.len() {
                ipsym = 0;
            }
        }

        let mut icon1 = vec![0i32; n];
        let mut icon2 = vec![0i32; n];

        for i in 0..n {
            let iz = (i + 1) as i32;
            let (xi1, yi1, zi1) = (x1[i], y1[i], z1[i]);
            let (xi2, yi2, zi2) = (x2[i], y2[i], z2[i]);
            let slen =
                ((xi2 - xi1) * (xi2 - xi1) + (yi2 - yi1) * (yi2 - yi1) + (zi2 - zi1) * (zi2 - zi1)).sqrt() * SMIN;

            // End 1.
            let mut on_ground = false;
            if ignd > 0 {
                if zi1 <= -slen {
                    return Err(Error::SegmentBelowGround(i + 1));
                }
                if zi1 <= slen {
                    icon1[i] = iz;
                    z1[i] = 0.0;
                    on_ground = true;
                }
            }
            if !on_ground {
                let mut ic = i;
                for _ in 1..n {
                    ic += 1;
                    if ic >= n {
                        ic = 0;
                    }
                    let sep = (xi1 - x1[ic]).abs() + (yi1 - y1[ic]).abs() + (zi1 - z1[ic]).abs();
                    if sep <= slen {
                        icon1[i] = -((ic + 1) as i32);
                        break;
                    }
                    let sep = (xi1 - x2[ic]).abs() + (yi1 - y2[ic]).abs() + (zi1 - z2[ic]).abs();
                    if sep <= slen {
                        icon1[i] = (ic + 1) as i32;
                        break;
                    }
                }
            }

            // End 2.
            if ignd > 0 {
                if zi2 <= -slen {
                    return Err(Error::SegmentBelowGround(i + 1));
                }
                if zi2 <= slen {
                    if icon1[i] == iz {
                        return Err(Error::SegmentInGroundPlane(i + 1));
                    }
                    icon2[i] = iz;
                    z2[i] = 0.0;
                    continue;
                }
            }
            let mut ic = i;
            for _ in 1..n {
                ic += 1;
                if ic >= n {
                    ic = 0;
                }
                let sep = (xi2 - x1[ic]).abs() + (yi2 - y1[ic]).abs() + (zi2 - z1[ic]).abs();
                if sep <= slen {
                    icon2[i] = (ic + 1) as i32;
                    break;
                }
                let sep = (xi2 - x2[ic]).abs() + (yi2 - y2[ic]).abs() + (zi2 - z2[ic]).abs();
                if sep <= slen {
                    icon2[i] = -((ic + 1) as i32);
                    break;
                }
            }
        }

        // A wire end on the centre of a patch: the patch is divided into
        // four around the end, at its place in the list. The four are
        // searched too. If one of them lies under the same wire end again,
        // the patch is too small for that wire: NEC-2 would go on dividing
        // without end. `made_for` holds the wire end a patch was made for.
        let mut made_for: Vec<Option<(usize, bool)>> = vec![None; patches.len()];
        let mut i = 0;
        while i < patches.len() {
            i += 1;
            let (xs, ys, zs) = (patches[i - 1].x, patches[i - 1].y, patches[i - 1].z);
            for iseg in 0..n {
                let slen =
                    ((x2[iseg] - x1[iseg]).abs() + (y2[iseg] - y1[iseg]).abs() + (z2[iseg] - z1[iseg]).abs()) * SMIN;
                let on_end1 = (x1[iseg] - xs).abs() + (y1[iseg] - ys).abs() + (z1[iseg] - zs).abs() <= slen;
                let on_end2 = (x2[iseg] - xs).abs() + (y2[iseg] - ys).abs() + (z2[iseg] - zs).abs() <= slen;
                if !(on_end1 || on_end2) {
                    continue;
                }
                let end = (iseg, !on_end1);
                if made_for[i - 1] == Some(end) || patches[i - 1].area <= 0.0 {
                    return Err(Error::PatchTooSmall { patch: i, segment: iseg + 1 });
                }
                if on_end1 {
                    icon1[iseg] = PCHCON + i as i32;
                } else {
                    icon2[iseg] = PCHCON + i as i32;
                }
                subph(&mut patches, &mut mp, i);
                made_for.splice(i - 1..i, [Some(end); 4]);
                break;
            }
        }

        if (n + patches.len()) / (np + mp) != 1 && ipsym == 0 {
            return Err(Error::Symmetry { cell: np + mp, total: n + patches.len() });
        }

        // Walking from a segment end through the connected ends must lead
        // back to that segment. The tolerance of a connection is relative to
        // the segment that is tested, so with very unequal segment lengths
        // the walk can enter a cycle that never returns; every later walk
        // over the junction (here and in the current expansion) relies on it.
        for j in 0..n {
            let jx = (j + 1) as i32;
            for end in [-1i32, 1] {
                let mut jend = end;
                let mut ix = if end < 0 { icon1[j] } else { icon2[j] };
                let mut steps = 0;
                // A patch ends the walk like a free end.
                while ix != 0 && ix.abs() != jx && ix <= PCHCON {
                    if ix < 0 {
                        ix = -ix;
                    } else {
                        jend = -jend;
                    }
                    let ixx = (ix - 1) as usize;
                    ix = if jend != 1 { icon1[ixx] } else { icon2[ixx] };
                    steps += 1;
                    if steps > 2 * n {
                        return Err(Error::JunctionNotClosed(j + 1));
                    }
                }
                if (ix == 0 || ix > PCHCON) && steps > 0 {
                    return Err(Error::JunctionNotClosed(j + 1));
                }
            }
        }

        // Make connected ends coincide exactly and list the junctions of
        // three or more segments. Each junction is handled once, from its
        // lowest-numbered segment.
        let mut junctions = Vec::new();
        let mut jco: Vec<i32> = Vec::new();
        for j in 0..n {
            let jx = (j + 1) as i32;
            // First end 1 of segment j, then end 2.
            for end in [-1i32, 1] {
                let mut jend = end;
                let mut ix = if end < 0 { icon1[j] } else { icon2[j] };
                jco.clear();
                jco.push(end * jx);
                let (mut xa, mut ya, mut za) = if end < 0 { (x1[j], y1[j], z1[j]) } else { (x2[j], y2[j], z2[j]) };

                if ix == 0 || ix == jx || ix > PCHCON {
                    continue;
                }
                let mut already_done = false;
                loop {
                    if ix < 0 {
                        ix = -ix;
                    } else {
                        jend = -jend;
                    }
                    if ix == jx {
                        break;
                    }
                    if ix < jx {
                        already_done = true;
                        break;
                    }
                    jco.push(ix * jend);
                    let ixx = (ix - 1) as usize;
                    if jend != 1 {
                        xa += x1[ixx];
                        ya += y1[ixx];
                        za += z1[ixx];
                        ix = icon1[ixx];
                    } else {
                        xa += x2[ixx];
                        ya += y2[ixx];
                        za += z2[ixx];
                        ix = icon2[ixx];
                    }
                    if ix == 0 {
                        break;
                    }
                }
                if already_done {
                    continue;
                }

                let count = jco.len() as f64;
                xa /= count;
                ya /= count;
                za /= count;
                for &entry in &jco {
                    if entry <= 0 {
                        let k = (-entry - 1) as usize;
                        x1[k] = xa;
                        y1[k] = ya;
                        z1[k] = za;
                    } else {
                        let k = (entry - 1) as usize;
                        x2[k] = xa;
                        y2[k] = ya;
                        z2[k] = za;
                    }
                }
                if jco.len() >= 3 {
                    junctions.push(jco.clone());
                }
            }
        }

        // Segment centres, lengths and direction cosines.
        let mut x = vec![0.0; n];
        let mut y = vec![0.0; n];
        let mut z = vec![0.0; n];
        let mut si = vec![0.0; n];
        let mut cab = vec![0.0; n];
        let mut sab = vec![0.0; n];
        let mut salp = vec![0.0; n];
        for i in 0..n {
            let xw1 = x2[i] - x1[i];
            let yw1 = y2[i] - y1[i];
            let zw1 = z2[i] - z1[i];
            x[i] = (x1[i] + x2[i]) / 2.0;
            y[i] = (y1[i] + y2[i]) / 2.0;
            z[i] = (z1[i] + z2[i]) / 2.0;
            let xw2 = xw1 * xw1 + yw1 * yw1 + zw1 * zw1;
            let mut yw2 = xw2.sqrt();
            yw2 = (xw2 / yw2 + yw2) * 0.5;
            si[i] = yw2;
            cab[i] = xw1 / yw2;
            sab[i] = yw1 / yw2;
            salp[i] = (zw1 / yw2).clamp(-1.0, 1.0);
            let valid = si[i].is_finite() && si[i] > 1.0e-20;
            if !valid {
                return Err(Error::SegmentData { segment: i + 1, tag: itag[i] });
            }
        }

        // Wires that occupy the same space.
        if overlap != OverlapCheck::Off {
            let ends = Ends { x1: &x1, y1: &y1, z1: &z1, x2: &x2, y2: &y2, z2: &z2 };
            let pairs = overlaps(&ends, [&x, &y, &z], &si, &bi, &itag);
            if !pairs.is_empty() {
                if overlap == OverlapCheck::Error {
                    return Err(Error::Overlap(pairs));
                }
                warnings.extend(overlap_lines(&pairs));
            }
        }

        Ok(Structure {
            np,
            mp,
            ipsym,
            ground_flag: ignd,
            patches,
            itag,
            icon1,
            icon2,
            x1,
            y1,
            z1,
            x2,
            y2,
            z2,
            x,
            y,
            z,
            si,
            bi,
            cab,
            sab,
            salp,
            junctions,
            warnings,
        })
    }
}

/// A finished structure: segments in metres with their connection data.
///
/// Arrays are indexed by segment (0-based). Connection numbers in `icon1` and
/// `icon2` follow NEC-2: 0 is a free end, `+k` / `-k` name the connected
/// segment `k` (1-based, the sign tells which of its ends is met) and the
/// segment's own number marks an end that connects to its image in the ground.
#[derive(Debug, Clone)]
pub struct Structure {
    /// Number of segments and of patches in one symmetric cell.
    pub np: usize,
    pub mp: usize,
    /// Symmetry flag: 0 none, positive for plane symmetry, negative for rotation.
    pub ipsym: i32,
    /// Ground-plane flag of the GE card.
    pub ground_flag: i32,
    /// Surface patches. A patch that a wire end is connected to has been
    /// divided into four.
    pub patches: Vec<Patch>,
    pub itag: Vec<i32>,
    pub icon1: Vec<i32>,
    pub icon2: Vec<i32>,
    pub x1: Vec<f64>,
    pub y1: Vec<f64>,
    pub z1: Vec<f64>,
    pub x2: Vec<f64>,
    pub y2: Vec<f64>,
    pub z2: Vec<f64>,
    /// Segment centres.
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    /// Segment lengths.
    pub si: Vec<f64>,
    /// Wire radii.
    pub bi: Vec<f64>,
    /// Direction cosines of the segments: cos(a)cos(b), cos(a)sin(b), sin(a).
    pub cab: Vec<f64>,
    pub sab: Vec<f64>,
    pub salp: Vec<f64>,
    /// Junctions of three or more segments: signed segment numbers, negative
    /// for end 1 and positive for end 2.
    pub junctions: Vec<Vec<i32>>,
    /// Remarks about the structure that did not stop the input.
    pub warnings: Vec<String>,
}

impl Structure {
    /// Number of segments.
    pub fn n(&self) -> usize {
        self.bi.len()
    }

    /// Number of patches.
    pub fn m(&self) -> usize {
        self.patches.len()
    }

    /// 1-based number of the `mx`-th segment with tag `itagi` (see [`isegno`]).
    pub fn isegno(&self, itagi: i32, mx: i32) -> Result<usize> {
        isegno(&self.itag, itagi, mx)
    }

    /// Number of symmetric cells the structure consists of.
    pub fn symmetry_sections(&self) -> usize {
        (self.n() + self.m()) / (self.np + self.mp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cone angle is the angle whose tangent is the change of the radius
    /// over the length, also where the formula of NEC-2 divides zero by zero.
    #[test]
    fn cone_angle_of_a_spiral() {
        let cone_angle = |a1: f64, a2: f64| {
            let info = Geometry::new().helix(1, 36, 0.12, 0.36, a1, a1, a2, a2, 0.001).unwrap();
            let Some(HelixInfo::Spiral { cone_angle }) = info else { panic!("{info:?} is no spiral") };
            cone_angle
        };
        for (a1, a2) in [(0.08, 0.02), (0.02, 0.08), (0.08, 0.0)] {
            let expected = f64::atan((a2 - a1) / 0.36);
            assert!((cone_angle(a1, a2) - expected).abs() < 1.0e-15, "{a1} to {a2}: {}", cone_angle(a1, a2));
        }
    }

    #[test]
    fn straight_wire_is_cut_into_equal_segments() {
        let mut g = Geometry::new();
        g.wire(1, 4, [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 0.001, 1.0, 1.0).unwrap();
        let s = g.connect(0).unwrap();
        assert_eq!(s.n(), 4);
        assert_eq!(s.icon1, vec![0, 1, 2, 3]);
        assert_eq!(s.icon2, vec![2, 3, 4, 0]);
        for i in 0..4 {
            assert!((s.si[i] - 0.5).abs() < 1e-15);
            assert_eq!(s.salp[i], 1.0);
        }
        assert_eq!(s.z, vec![-0.75, -0.25, 0.25, 0.75]);
    }

    #[test]
    fn wires_that_start_at_the_same_point_are_joined_end1_to_end1() {
        let mut g = Geometry::new();
        g.wire(1, 1, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.001, 1.0, 1.0).unwrap();
        g.wire(2, 1, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.001, 1.0, 1.0).unwrap();
        let s = g.connect(0).unwrap();
        assert_eq!(s.icon1, vec![-2, -1]);
        assert_eq!(s.icon2, vec![0, 0]);
        assert!(s.junctions.is_empty());
    }

    #[test]
    fn wire_base_on_ground_connects_to_its_image() {
        let mut g = Geometry::new();
        g.wire(1, 2, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.001, 1.0, 1.0).unwrap();
        let s = g.connect(1).unwrap();
        assert_eq!(s.icon1, vec![1, 1]);
        assert_eq!(s.icon2, vec![2, 0]);
    }

    #[test]
    fn patches_of_the_four_shapes() {
        let mut g = Geometry::new();
        // Any shape: centre, normal straight up, area.
        g.patch(0, 1, [1.0, 2.0, 3.0], [PI / 2.0, 0.0, 0.5], [0.0; 3], [0.0; 3]).unwrap();
        // Rectangle, triangle and quadrilateral in the plane z = 0.
        let corners = ([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        g.patch(0, 2, corners.0, corners.1, corners.2, corners.3).unwrap();
        g.patch(0, 3, corners.0, corners.1, corners.2, corners.3).unwrap();
        g.patch(0, 4, corners.0, corners.1, corners.2, corners.3).unwrap();
        // A rectangle of 2 by 3 patches.
        g.patch(2, 3, corners.0, corners.1, [2.0, 3.0, 0.0], [0.0; 3]).unwrap();
        // Corners that do not lie in a plane.
        assert_eq!(g.patch(0, 4, corners.0, corners.1, corners.2, [0.0, 1.0, 1.0]), Err(Error::PatchNotPlanar));

        let s = g.connect(0).unwrap();
        assert_eq!(s.m(), 4 + 6);
        let areas: Vec<f64> = s.patches.iter().map(|p| p.area).collect();
        assert_eq!(areas, [0.5, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        assert_eq!((s.patches[0].x, s.patches[0].y, s.patches[0].z), (1.0, 2.0, 3.0));
        assert_eq!((s.patches[1].x, s.patches[1].y), (1.0, 0.5));
        assert!((s.patches[2].x - 4.0 / 3.0).abs() < 1e-15 && (s.patches[2].y - 1.0 / 3.0).abs() < 1e-15);
        assert_eq!((s.patches[9].x, s.patches[9].y), (1.5, 2.5));
        for patch in &s.patches {
            let [nx, ny, nz] = patch.normal();
            assert!(nx.abs() < 1e-8 && ny.abs() < 1e-8 && (nz - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn a_wire_end_divides_its_patch_into_four() {
        let mut g = Geometry::new();
        g.patch(3, 1, [-0.3, -0.1, 0.0], [0.3, -0.1, 0.0], [0.3, 0.1, 0.0], [0.0; 3]).unwrap();
        g.wire(1, 2, [0.0, 0.0, 0.0], [0.0, 0.0, 0.2], 0.001, 1.0, 1.0).unwrap();
        let s = g.connect(0).unwrap();
        // The middle patch is replaced by four quarters at its place.
        assert_eq!((s.m(), s.mp), (6, 6));
        assert_eq!(s.icon1[0], PCHCON + 2);
        let signs = [(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)];
        for (quarter, (sx, sy)) in s.patches[1..5].iter().zip(signs) {
            assert!((quarter.x - 0.05 * sx).abs() < 1e-15 && (quarter.y - 0.05 * sy).abs() < 1e-15);
            assert!((quarter.area - 0.01).abs() < 1e-15);
        }
        assert!((s.patches[0].x + 0.2).abs() < 1e-15 && (s.patches[5].x - 0.2).abs() < 1e-15);
    }

    #[test]
    fn segment_selection_by_tag() {
        let tags = [1, 1, 2, 2, 2];
        assert_eq!(isegno(&tags, 2, 2), Ok(4));
        assert_eq!(isegno(&tags, 0, 5), Ok(5));
        assert_eq!(isegno(&tags, 3, 1), Err(Error::NoSegmentWithTag(3)));
        assert_eq!(isegno(&tags, 1, 0), Err(Error::ZeroSegmentPosition));
        assert_eq!(isegno(&tags, 0, 6), Err(Error::SegmentOutOfRange(6)));
    }
}
