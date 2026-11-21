//! Radiated fields: far-field patterns, gains and the ground wave.
//!
//! Ported from `ffld`, `gfld` and `rdpat` (`radiation.c`) of nec2c 1.3.1.

use num_complex::Complex64;

use super::{C0, Engine};
use crate::consts::{CONST3, CONST4, ETA, PI, TA, TD, TP, cang, cmplx, db10};

/// Parameters of a radiation-pattern calculation (RP card).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PatternRequest {
    /// Mode: 0 space wave, 1 space plus ground wave, 2 linear cliff,
    /// 3 circular cliff, 4 radial-wire screen, 5 screen and linear cliff,
    /// 6 screen and circular cliff.
    pub ifar: i32,
    /// Number of theta (or z) and phi values.
    pub nth: i32,
    pub nph: i32,
    /// Directive gain (true) instead of power gain.
    pub directive: bool,
    /// Average gain: 0 off, 1 computed, 2 computed without listing the pattern.
    pub iavp: i32,
    /// Normalized gain: 0 off, 1 major axis, 2 minor axis, 3 vertical,
    /// 4 horizontal, 5 total.
    pub inor: i32,
    /// First angles and increments, degrees (for mode 1: heights in metres).
    pub thets: f64,
    pub phis: f64,
    pub dth: f64,
    pub dph: f64,
    /// Range in metres; zero leaves the factor exp(-jkR)/R out.
    pub rfld: f64,
    /// Gain normalization factor in dB; zero takes the maximum gain.
    pub gnor: f64,
    /// Excitation type of the solution (EX card).
    pub ixtyp: i32,
    /// Sixth floating-point field of the EX card (axial ratio or current moment).
    pub xpr6: f64,
    /// Second medium: relative permittivity, conductivity, distance of the
    /// cliff edge and its height in metres (GN or GD card).
    pub epsr2: f64,
    pub sig2: f64,
    pub clt: f64,
    pub cht: f64,
}

/// One line of a far-field pattern.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FarFieldRow {
    pub theta: f64,
    pub phi: f64,
    /// Gains in dB of the major and minor polarization axes, of the vertical
    /// (theta) and horizontal (phi) components, and of the total field.
    pub gain_major: f64,
    pub gain_minor: f64,
    pub gain_vertical: f64,
    pub gain_horizontal: f64,
    pub gain_total: f64,
    pub axial_ratio: f64,
    /// Tilt of the polarization ellipse, degrees.
    pub tilt: f64,
    /// Sense of rotation: `""` (no field), `"LINEAR"`, `"RIGHT "` or `"LEFT  "`.
    pub sense: &'static str,
    /// Field magnitudes in volts per metre (at the requested range, or
    /// normalized to 1 m) and phases in degrees.
    pub e_theta_mag: f64,
    pub e_theta_phase: f64,
    pub e_phi_mag: f64,
    pub e_phi_phase: f64,
}

/// One line of the fields near the ground (RP mode 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundWaveRow {
    /// Cylindrical coordinates of the field point: metres and degrees.
    pub rho: f64,
    pub phi: f64,
    pub z: f64,
    pub e_theta_mag: f64,
    pub e_theta_phase: f64,
    pub e_phi_mag: f64,
    pub e_phi_phase: f64,
    pub e_radial_mag: f64,
    pub e_radial_phase: f64,
}

/// One entry of the normalized-gain table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GainSample {
    /// Direction in degrees.
    pub theta: f64,
    pub phi: f64,
    /// Gain in dB before normalization.
    pub gain: f64,
}

/// Gains for the normalized-gain table.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedGain {
    /// Which gain was stored (the `inor` of the request).
    pub kind: i32,
    /// Normalization factor in dB.
    pub factor: f64,
    /// The directions that were computed, theta running fastest. With a
    /// ground the directions below it are not among them.
    pub samples: Vec<GainSample>,
}

/// Average power gain over the directions of a pattern.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AverageGain {
    /// The power gain averaged over the solid angle that the samples cover,
    /// or `None` if they cover none.
    pub gain: Option<f64>,
    /// Solid angle that the samples cover, steradians.
    pub solid_angle: f64,
}

impl AverageGain {
    /// The part of the input power that is radiated into the sampled solid
    /// angle: the average gain times the sampled fraction of the sphere. For
    /// a pattern that covers all directions of radiation this estimates the
    /// efficiency of the antenna, and a value far from the expected one
    /// shows that the model or the sampling is not sound.
    pub fn efficiency(&self) -> Option<f64> {
        self.gain.map(|gain| gain * self.solid_angle / (4.0 * PI))
    }
}

/// Result of [`Engine::radiation_pattern`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RadiationPattern {
    pub rows: Vec<FarFieldRow>,
    pub ground_wave: Vec<GroundWaveRow>,
    /// `1 / R` and the phase of exp(-jkR) in degrees, when a range was given.
    pub range_factor: Option<(f64, f64)>,
    pub average_gain: Option<AverageGain>,
    pub normalized: Option<NormalizedGain>,
}

/// Integral of |sin(theta)| between two angles in radians, in any order.
fn sine_weight(a: f64, b: f64) -> f64 {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let (klo, khi) = ((lo / PI).floor(), (hi / PI).floor());
    if klo == khi {
        // The sine keeps its sign: the weight NEC-2 uses.
        return (lo.cos() - hi.cos()).abs();
    }
    // Up to the first zero of the sine, the whole lobes between the zeros,
    // and from the last zero on.
    (lo.cos() - ((klo + 1.0) * PI).cos()).abs() + 2.0 * (khi - klo - 1.0) + ((khi * PI).cos() - hi.cos()).abs()
}

/// sin(x)/x with the series NEC-2 uses for small arguments.
#[inline]
fn sinc(x: f64) -> f64 {
    if x.abs() >= 1.0e-7 { x.sin() / x } else { 1.0 - x * x / 6.0 }
}

impl Engine {
    /// Integrals of the three current terms of segment `i` against the phase
    /// factor for direction cosine `omega` along the segment: returns the
    /// real and imaginary parts of the segment's current moment.
    #[inline]
    fn current_moment(&self, i: usize, omega: f64) -> (f64, f64) {
        let el = PI * self.si[i];
        let sill = omega * el;
        let top = el + sill;
        let bot = el - sill;
        let a =
            if omega.abs() >= 1.0e-7 { 2.0 * sill.sin() / omega } else { (2.0 - omega * omega * el * el / 3.0) * el };
        let too = sinc(top);
        let boo = sinc(bot);
        let b = el * (boo - too);
        let c = el * (boo + too);
        let rr = a * self.air[i] + b * self.bii[i] + c * self.cir[i];
        let ri = a * self.aii[i] - b * self.bir[i] + c * self.cii[i];
        (rr, ri)
    }

    /// Far-zone radiated field in direction (`thet`, `phi`), radians; the
    /// factor exp(-jkr)/(r/lambda) is not included. Returns (E theta, E phi).
    pub(crate) fn ffld(&self, thet: f64, phi: f64) -> (Complex64, Complex64) {
        let gnd = &self.gnd;
        let phx = -phi.sin();
        let phy = phi.cos();
        let mut roz = thet.cos();
        let thx = roz * phy;
        let thy = -roz * phx;
        let thz = -thet.sin();
        let rox = -thz * phy;
        let roy = thz * phx;

        let (mut cix, mut ciy, mut ciz) = (C0, C0, C0);
        let (mut ccx, mut ccy, mut ccz) = (C0, C0, C0);
        let (mut rrv, mut rrh) = (C0, C0);
        let (mut rrv1, mut rrh1, mut rrv2, mut rrh2) = (C0, C0, C0, C0);
        let (mut tthet, mut darg) = (0.0, 0.0);

        // The structure, then its image in the ground.
        for k in 0..gnd.ksymp {
            if k != 0 {
                if gnd.iperf == 1 {
                    rrv = cmplx(-1.0, 0.0);
                    rrh = cmplx(-1.0, 0.0);
                } else {
                    // Fresnel reflection coefficients of the ground.
                    let zrsin = (1.0 - gnd.zrati * gnd.zrati * thz * thz).sqrt();
                    rrv = -(roz - gnd.zrati * zrsin) / (roz + gnd.zrati * zrsin);
                    rrh = (gnd.zrati * roz - zrsin) / (gnd.zrati * roz + zrsin);
                }

                // The cliff problem needs the coefficients of both media.
                if gnd.ifar > 1 {
                    rrv1 = rrv;
                    rrh1 = rrh;
                    tthet = thet.tan();
                    if gnd.ifar != 4 {
                        let zrsin = (1.0 - gnd.zrati2 * gnd.zrati2 * thz * thz).sqrt();
                        rrv2 = -(roz - gnd.zrati2 * zrsin) / (roz + gnd.zrati2 * zrsin);
                        rrh2 = (gnd.zrati2 * roz - zrsin) / (gnd.zrati2 * roz + zrsin);
                        darg = -TP * 2.0 * gnd.ch * roz;
                    }
                }

                roz = -roz;
                ccx = cix;
                ccy = ciy;
                ccz = ciz;
            }

            cix = C0;
            ciy = C0;
            ciz = C0;

            for i in 0..self.n {
                let omega = -(rox * self.cab[i] + roy * self.sab[i] + roz * self.salp[i]);
                let (rr, ri) = self.current_moment(i, omega);
                let mut arg = TP * (self.x[i] * rox + self.y[i] * roy + self.z[i] * roz);

                if k != 1 || gnd.ifar < 2 {
                    let exa = cmplx(arg.cos(), arg.sin()) * cmplx(rr, ri);
                    cix += exa * self.cab[i];
                    ciy += exa * self.sab[i];
                    ciz += exa * self.salp[i];
                    continue;
                }

                // Image contribution in the cliff and ground-screen problems:
                // the reflection coefficient depends on where the ray is
                // reflected. `dr` is the distance of the specular point.
                let dr = self.z[i] * tthet;
                let mut d = dr * phy + self.x[i];
                let second_medium;
                if gnd.ifar == 2 {
                    second_medium = gnd.cl - d <= 0.0;
                } else {
                    d = (d * d + (self.y[i] - dr * phx) * (self.y[i] - dr * phx)).sqrt();
                    if gnd.ifar == 3 {
                        second_medium = gnd.cl - d <= 0.0;
                    } else if gnd.nradl > 0 && gnd.scrwl - d >= 0.0 {
                        // Reflection from the radial-wire ground screen.
                        d += gnd.t2;
                        let mut zscrn = gnd.t1 * d * (d / gnd.t2).ln();
                        zscrn = (zscrn * gnd.zrati) / (ETA * gnd.zrati + zscrn);
                        let zrsin = (1.0 - zscrn * zscrn * thz * thz).sqrt();
                        rrv = (roz + zscrn * zrsin) / (-roz + zscrn * zrsin);
                        rrh = (zscrn * roz + zrsin) / (zscrn * roz - zrsin);
                        let exa = cmplx(arg.cos(), arg.sin()) * cmplx(rr, ri);
                        self.add_image(exa, i, phx, phy, rrv, rrh, &mut cix, &mut ciy, &mut ciz);
                        continue;
                    } else if gnd.ifar == 4 {
                        second_medium = false;
                    } else {
                        if gnd.ifar == 5 {
                            d = dr * phy + self.x[i];
                        }
                        second_medium = gnd.cl - d <= 0.0;
                    }
                }

                if second_medium {
                    rrv = rrv2;
                    rrh = rrh2;
                    arg += darg;
                } else {
                    rrv = rrv1;
                    rrh = rrh1;
                }
                let exa = cmplx(arg.cos(), arg.sin()) * cmplx(rr, ri);
                self.add_image(exa, i, phx, phy, rrv, rrh, &mut cix, &mut ciy, &mut ciz);
            }

            if k == 0 {
                continue;
            }

            // Image of the structure in an infinite planar ground.
            if gnd.ifar < 2 {
                let cdp = (cix * phx + ciy * phy) * (rrh - rrv);
                cix = ccx + cix * rrv + cdp * phx;
                ciy = ccy + ciy * rrv + cdp * phy;
                ciz = ccz - ciz * rrv;
            } else {
                cix += ccx;
                ciy += ccy;
                ciz += ccz;
            }
        }

        if self.m == 0 {
            let eth = (cix * thx + ciy * thy + ciz * thz) * CONST3;
            let eph = (cix * phx + ciy * phy) * CONST3;
            return (eth, eph);
        }

        // The field of the patch currents, then of their image in the ground.
        let roz = thet.cos();
        let mut rfl = -1.0;
        let (mut ex, mut ey, mut ez) = (C0, C0, C0);
        for ip in 0..gnd.ksymp {
            rfl = -rfl;
            let [mut gx, mut gy, mut gz] = self.fflds(rox, roy, roz * rfl);
            if ip != 1 {
                ex = gx;
                ey = gy;
                ez = gz;
                continue;
            }
            if gnd.iperf == 1 {
                gx = -gx;
                gy = -gy;
                gz = -gz;
            } else {
                let mut rrv = (1.0 - gnd.zrati * gnd.zrati * thz * thz).sqrt();
                let mut rrh = gnd.zrati * roz;
                rrh = (rrh - rrv) / (rrh + rrv);
                rrv = gnd.zrati * rrv;
                rrv = -(roz - rrv) / (roz + rrv);
                let across = (gx * phx + gy * phy) * (rrh - rrv);
                gx = gx * rrv + across * phx;
                gy = gy * rrv + across * phy;
                gz *= rrv;
            }
            ex += gx;
            ey += gy;
            ez -= gz;
        }
        ex += cix * CONST3;
        ey += ciy * CONST3;
        ez += ciz * CONST3;
        (ex * thx + ey * thy + ez * thz, ex * phx + ey * phy)
    }

    /// The x, y and z components of the far electric field of the patch
    /// currents in the direction `(rox, roy, roz)` (`fflds`).
    fn fflds(&self, rox: f64, roy: f64, roz: f64) -> [Complex64; 3] {
        let scur = &self.cur[self.n..];
        let (mut ex, mut ey, mut ez) = (C0, C0, C0);
        for (j, p) in self.patches.iter().enumerate() {
            let arg = TP * (rox * p.x + roy * p.y + roz * p.z);
            let ct = cmplx(arg.cos() * p.area, arg.sin() * p.area);
            ex += scur[3 * j] * ct;
            ey += scur[3 * j + 1] * ct;
            ez += scur[3 * j + 2] * ct;
        }
        let ct = rox * ex + roy * ey + roz * ez;
        [CONST4 * (ct * rox - ex), CONST4 * (ct * roy - ey), CONST4 * (ct * roz - ez)]
    }

    /// Adds the contribution of the image of segment `i`, modified by the
    /// reflection coefficients `rrv` and `rrh`.
    #[allow(clippy::too_many_arguments)]
    #[inline]
    fn add_image(
        &self,
        exa: Complex64,
        i: usize,
        phx: f64,
        phy: f64,
        rrv: Complex64,
        rrh: Complex64,
        cix: &mut Complex64,
        ciy: &mut Complex64,
        ciz: &mut Complex64,
    ) {
        let tix = exa * self.cab[i];
        let tiy = exa * self.sab[i];
        let tiz = exa * self.salp[i];
        let cdp = (tix * phx + tiy * phy) * (rrh - rrv);
        *cix = *cix + tix * rrv + cdp * phx;
        *ciy = *ciy + tiy * rrv + cdp * phy;
        *ciz -= tiz * rrv;
    }

    /// Radiated field including the ground wave at the point with cylindrical
    /// coordinates `rho`, `phi`, `rz` (wavelengths, radians). Returns
    /// (E theta, E phi, E radial).
    pub(crate) fn gfld(&mut self, rho: f64, phi: f64, rz: f64, ux: Complex64) -> (Complex64, Complex64, Complex64) {
        let r = (rho * rho + rz * rz).sqrt();
        if self.gnd.ksymp == 1 || ux.norm() > 0.5 || r > 1.0e5 {
            // Space wave only.
            let thet = if rz >= 1.0e-20 { (rho / rz).atan() } else { PI * 0.5 };
            let (eth, epi) = self.ffld(thet, phi);
            let arg = -TP * r;
            let exa = cmplx(arg.cos(), arg.sin()) / r;
            return (eth * exa, epi * exa, C0);
        }

        // Space and ground waves.
        self.gwav.u = ux;
        self.gwav.u2 = ux * ux;
        let phx = -phi.sin();
        let phy = phi.cos();
        let rx = rho * phy;
        let ry = -rho * phx;
        let (mut cix, mut ciy, mut ciz) = (C0, C0, C0);

        for i in 0..self.n {
            let dx = self.cab[i];
            let dy = self.sab[i];
            let dz = self.salp[i];
            let rix = rx - self.x[i];
            let riy = ry - self.y[i];
            let rhs = rix * rix + riy * riy;
            let rhp = rhs.sqrt();
            let (rhx, rhy) = if rhp >= 1.0e-6 { (rix / rhp, riy / rhp) } else { (1.0, 0.0) };

            let mut calp = 1.0 - dz * dz;
            let (cph, sph);
            if calp >= 1.0e-6 {
                calp = calp.sqrt();
                let cbet = dx / calp;
                let sbet = dy / calp;
                cph = rhx * cbet + rhy * sbet;
                sph = rhy * cbet - rhx * sbet;
            } else {
                cph = rhx;
                sph = rhy;
            }

            // Current moment of the segment and of its image.
            let mut rfl = -1.0;
            for k in 0..2 {
                rfl = -rfl;
                let riz = rz - self.z[i] * rfl;
                let rxyz = (rix * rix + riy * riy + riz * riz).sqrt();
                let rnx = rix / rxyz;
                let rny = riy / rxyz;
                let rnz = riz / rxyz;
                let omega = -(rnx * dx + rny * dy + rnz * dz * rfl);
                let (rr, ri) = self.current_moment(i, omega);
                let arg = TP * (self.x[i] * rnx + self.y[i] * rny + self.z[i] * rnz * rfl);
                let exa = cmplx(arg.cos(), arg.sin()) * cmplx(rr, ri) / TP;
                if k != 1 {
                    self.gwav.xx1 = exa;
                    self.gwav.r1 = rxyz;
                    self.gwav.zmh = riz;
                } else {
                    self.gwav.xx2 = exa;
                    self.gwav.r2 = rxyz;
                    self.gwav.zph = riz;
                }
            }

            let (erv, ezv, mut erh, mut ezh, mut eph) = self.gwave();
            erh = erh * cph * calp + erv * dz;
            eph = eph * sph * calp;
            ezh = ezh * cph * calp + ezv * dz;
            let ex = erh * rhx - eph * rhy;
            let ey = erh * rhy + eph * rhx;
            cix += ex;
            ciy += ey;
            ciz += ezh;
        }

        let arg = -TP * r;
        let exa = cmplx(arg.cos(), arg.sin());
        cix *= exa;
        ciy *= exa;
        ciz *= exa;
        let rnx = rx / r;
        let rny = ry / r;
        let rnz = rz / r;
        let thx = rnz * phy;
        let thy = -rnz * phx;
        let thz = -rho / r;
        let eth = cix * thx + ciy * thy + ciz * thz;
        let epi = cix * phx + ciy * phy;
        let erd = cix * rnx + ciy * rny + ciz * rnz;
        (eth, epi, erd)
    }

    /// Computes the radiation pattern, gains and normalized gain for the
    /// last solution (`rdpat`).
    ///
    /// With a ground, a direction below it is not a sample: it is left out
    /// of the pattern, of the normalized gains and of the average gain.
    pub fn radiation_pattern(&mut self, req: &PatternRequest) -> RadiationPattern {
        let mut out = RadiationPattern::default();
        let wlam = self.wlam;
        self.gnd.ifar = req.ifar;
        let ifar = req.ifar;
        let nth = req.nth.max(0) as usize;
        let mut samples = Vec::new();

        // The theta values (heights in mode 1), accumulated as NEC-2 does.
        let mut thetas = Vec::with_capacity(nth);
        let mut thet = req.thets - req.dth;
        for _ in 0..nth {
            thet += req.dth;
            thetas.push(thet);
        }
        // NEC-2 leaves out theta above 90.01 degrees over a ground. The test
        // on the cosine does the same for angles of any sign and size.
        let ground = self.gnd.ksymp == 2 && ifar != 1;
        let horizon = (90.01 * TA).cos();
        let valid: Vec<bool> = thetas.iter().map(|&thet| !(ground && (thet * TA).cos() < horizon)).collect();

        // Second medium for the cliff problem.
        if ifar > 1 && ifar != 4 {
            self.gnd.cl = req.clt / wlam;
            self.gnd.ch = req.cht / wlam;
            self.gnd.zrati2 = (1.0 / cmplx(req.epsr2, -req.sig2 * wlam * 59.96)).sqrt();
        }

        let (mut exrm, mut exra) = (0.0, 0.0);
        if ifar != 1 && req.rfld >= 1.0e-20 {
            exrm = 1.0 / req.rfld;
            exra = req.rfld / wlam;
            exra = -360.0 * (exra - exra.floor());
            out.range_factor = Some((exrm, exra));
        }

        // Gain normalization.
        let pinr = if req.ixtyp == 4 { 394.51 * req.xpr6 * req.xpr6 * wlam * wlam } else { self.pin };
        let (gcon, gcop);
        if req.ixtyp == 0 || req.ixtyp == 5 || req.ixtyp == 4 {
            gcop = wlam * wlam * 2.0 * PI / (376.73 * pinr);
            let prad = pinr - self.structure_loss() - self.pnls;
            gcon = if req.directive { gcop * pinr / prad } else { gcop };
        } else {
            gcon = 4.0 * PI / (1.0 + req.xpr6 * req.xpr6);
            gcop = gcon;
        }

        let mut gmax = -1.0e+10;
        let mut pint = 0.0;
        let mut solid = 0.0;
        let tmp1 = req.dph * TA;
        let tmp2 = 0.5 * req.dth * TA;
        let mut phi = req.phis - req.dph;

        for kph in 1..=req.nph {
            phi += req.dph;
            let pha = phi * TA;

            for (kth, &thet) in thetas.iter().enumerate() {
                if !valid[kth] {
                    continue;
                }
                let tha = thet * TA;

                if ifar == 1 {
                    let zrati = self.gnd.zrati;
                    let (eth, eph, erd) = self.gfld(req.rfld / wlam, pha, thet / wlam, zrati);
                    out.ground_wave.push(GroundWaveRow {
                        rho: req.rfld,
                        phi,
                        z: thet,
                        e_theta_mag: eth.norm_sqr().sqrt(),
                        e_theta_phase: cang(eth),
                        e_phi_mag: eph.norm_sqr().sqrt(),
                        e_phi_phase: cang(eph),
                        e_radial_mag: erd.norm(),
                        e_radial_phase: cang(erd),
                    });
                    continue;
                }

                let (eth, eph) = self.ffld(tha, pha);
                let ethm2 = eth.norm_sqr();
                let mut ethm = ethm2.sqrt();
                let mut etha = cang(eth);
                let ephm2 = eph.norm_sqr();
                let mut ephm = ephm2.sqrt();
                let mut epha = cang(eph);

                // Polarization ellipse.
                let (tilta, emajr2, eminr2, axrat, isens);
                if ethm2 <= 1.0e-20 && ephm2 <= 1.0e-20 {
                    tilta = 0.0;
                    emajr2 = 0.0;
                    eminr2 = 0.0;
                    axrat = 0.0;
                    isens = "";
                } else {
                    let mut dfaz = epha - etha;
                    let dfaz2 = if epha >= 0.0 { dfaz - 360.0 } else { dfaz + 360.0 };
                    if dfaz.abs() > dfaz2.abs() {
                        dfaz = dfaz2;
                    }
                    let cdfaz = (dfaz * TA).cos();
                    let mut tstor1 = ethm2 - ephm2;
                    let mut tstor2 = 2.0 * ephm * ethm * cdfaz;
                    let tilt = 0.5 * tstor2.atan2(tstor1);
                    let stilta = tilt.sin();
                    tstor1 = tstor1 * stilta * stilta;
                    tstor2 = tstor2 * stilta * tilt.cos();
                    emajr2 = -tstor1 + tstor2 + ethm2;
                    eminr2 = (tstor1 - tstor2 + ephm2).max(0.0);
                    axrat = (eminr2 / emajr2).sqrt();
                    tilta = tilt * TD;
                    isens = if axrat <= 1.0e-5 {
                        "LINEAR"
                    } else if dfaz <= 0.0 {
                        "RIGHT "
                    } else {
                        "LEFT  "
                    };
                }

                let gnmj = db10(gcon * emajr2);
                let gnmn = db10(gcon * eminr2);
                let gnv = db10(gcon * ethm2);
                let gnh = db10(gcon * ephm2);
                let gtot = db10(gcon * (ethm2 + ephm2));

                if req.inor > 0 {
                    let stored = match req.inor {
                        1 => gnmj,
                        2 => gnmn,
                        3 => gnv,
                        4 => gnh,
                        _ => gtot,
                    };
                    // The angles of the table are not accumulated in NEC-2.
                    samples.push(GainSample {
                        theta: req.thets + kth as f64 * req.dth,
                        phi: req.phis + f64::from(kph - 1) * req.dph,
                        gain: stored,
                    });
                    if stored > gmax {
                        gmax = stored;
                    }
                }

                if req.iavp != 0 {
                    let power = gcop * (ethm2 + ephm2);
                    // Solid angle of the sample: its cell reaches half a
                    // step towards each neighbour that is a sample too.
                    let tmp3 = if kth > 0 && valid[kth - 1] { tha - tmp2 } else { tha };
                    let tmp4 = if kth + 1 < nth && valid[kth + 1] { tha + tmp2 } else { tha };
                    let mut da = (tmp1 * sine_weight(tmp3, tmp4)).abs();
                    if kph == 1 || kph == req.nph {
                        da *= 0.5;
                    }
                    pint += power * da;
                    solid += da;
                    if req.iavp == 2 {
                        continue;
                    }
                }

                ethm *= wlam;
                ephm *= wlam;
                if req.rfld >= 1.0e-20 {
                    ethm *= exrm;
                    etha += exra;
                    ephm *= exrm;
                    epha += exra;
                }

                out.rows.push(FarFieldRow {
                    theta: thet,
                    phi,
                    gain_major: gnmj,
                    gain_minor: gnmn,
                    gain_vertical: gnv,
                    gain_horizontal: gnh,
                    gain_total: gtot,
                    axial_ratio: axrat,
                    tilt: tilta,
                    sense: isens,
                    e_theta_mag: ethm,
                    e_theta_phase: etha,
                    e_phi_mag: ephm,
                    e_phi_phase: epha,
                });
            }
        }

        if req.iavp != 0 {
            // The gain is averaged over the solid angle that was sampled
            // (NEC-2 divides by the nominal solid angle of the angle ranges).
            let average = AverageGain { gain: (solid >= 1.0e-20).then(|| pint / solid), solid_angle: solid };
            // An antenna cannot radiate more than its input power.
            let antenna = matches!(req.ixtyp, 0 | 4 | 5);
            if antenna && average.efficiency().is_some_and(|fraction| fraction > 1.01) {
                self.diag.warnings.push(
                    "WARNING: THE AVERAGE POWER GAIN EXCEEDS ITS LIMIT OF 4*PI DIVIDED BY THE SOLID ANGLE BY MORE \
                     THAN 1 PERCENT - THE PATTERN IS SAMPLED TOO COARSELY OR THE MODEL IS NOT VALID"
                        .to_string(),
                );
            }
            out.average_gain = Some(average);
        }

        // Without a sample there is nothing to normalize.
        if req.inor > 0 && !samples.is_empty() {
            if req.gnor.abs() > 1.0e-20 {
                gmax = req.gnor;
            }
            out.normalized = Some(NormalizedGain { kind: req.inor, factor: gmax, samples });
        }
        out
    }
}
