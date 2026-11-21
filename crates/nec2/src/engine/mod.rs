//! The NEC-2 solution engine.
//!
//! [`Engine`] holds the state that NEC-2 keeps in its COMMON blocks (and nec2c
//! in global structs): the structure scaled to wavelengths, the interaction
//! matrix, the current expansion and the ground parameters. Its methods are the
//! NEC-2 routines, kept under their original names (`tbf`, `efld`, `cmset`,
//! `netwk`, `rdpat` and so on) so that the code can be read next to the NEC-2
//! documentation (Part II, program description).
//!
//! The order of calls for one solution is the order of NEC-2 itself:
//!
//! 1. [`Engine::set_frequency`] scales the structure to wavelengths;
//! 2. [`Engine::set_loads`] and [`Engine::set_ground`] describe the environment;
//! 3. [`Engine::fill_and_factor`] builds and factors the interaction matrix;
//! 4. [`Engine::excite`] fills the right-hand side;
//! 5. [`Engine::solve`] solves for the currents, with networks if any;
//! 6. [`Engine::radiation_pattern`] and the other output routines use the currents.

use num_complex::Complex64;

use crate::error::{Error, Result};
use crate::geometry::{Patch, Structure};

mod basis;
mod excitation;
mod ground;
mod kernel;
mod loading;
mod matrix;
mod nearfield;
mod network;
mod patch;
mod radiation;
mod somnec;

pub use excitation::Excitation;
pub use ground::{GroundInfo, GroundParams, GroundType};
pub use loading::{Load, LoadLine, LoadingReport};
pub use network::{Asymmetry, Network, NetworkKind, PortRecord, SolveReport};
pub use radiation::{
    AverageGain, FarFieldRow, GainSample, GroundWaveRow, NormalizedGain, PatternRequest, RadiationPattern,
};

pub(crate) const C0: Complex64 = Complex64::new(0.0, 0.0);
pub(crate) const C1: Complex64 = Complex64::new(1.0, 0.0);

/// A smallest-to-largest pivot ratio below this value is reported as an
/// ill-conditioned matrix.
const ILL_CONDITIONED: f64 = 1.0e-5;

/// What the numerical routines report while they work.
#[derive(Debug, Clone, Default)]
pub(crate) struct Diag {
    /// Messages for the user, in the order they arose.
    pub warnings: Vec<String>,
    /// The adaptive integration that met a value that is not finite, or that
    /// did not end within its bound on the number of steps.
    pub fault: Option<&'static str>,
}

/// Bound on the passes of an adaptive integration that may halve its step
/// down to 1/`nm` of the interval. An integration that advances needs at
/// most `nm` steps and about as many changes of the step, so reaching the
/// bound means that the integration variable is not a number.
pub(crate) const fn pass_bound(nm: i32) -> u32 {
    4 * nm as u32
}

/// Whether all the values are finite.
pub(crate) fn all_finite(values: &[Complex64]) -> bool {
    values.iter().all(|v| v.re.is_finite() && v.im.is_finite())
}

/// The source segment of the field routines (COMMON /DATAJ/).
#[derive(Debug, Clone, Default)]
pub(crate) struct SourceSeg {
    /// Extended thin-wire kernel requested.
    pub iexk: i32,
    /// Extended-kernel status of end 1 and end 2: 0 use it, 1 free end, 2 do not use it.
    pub ind1: i32,
    pub ind2: i32,
    pub s: f64,
    pub b: f64,
    pub xj: f64,
    pub yj: f64,
    pub zj: f64,
    pub cabj: f64,
    pub sabj: f64,
    pub salpj: f64,
    /// For a source patch: its tangent vectors (the centre is in `xj`, `yj`,
    /// `zj` and the area in `s`), and 1 for the patch itself or 2 for its
    /// image in the ground.
    pub t1xj: f64,
    pub t1yj: f64,
    pub t1zj: f64,
    pub t2xj: f64,
    pub t2yj: f64,
    pub t2zj: f64,
    pub ipgnd: i32,
    /// Separation, in wavelengths, beyond which a segment is treated as a point source.
    pub rkh: f64,
    /// Field of the constant (k), sine (s) and cosine (c) current terms.
    pub exk: Complex64,
    pub eyk: Complex64,
    pub ezk: Complex64,
    pub exs: Complex64,
    pub eys: Complex64,
    pub ezs: Complex64,
    pub exc: Complex64,
    pub eyc: Complex64,
    pub ezc: Complex64,
}

/// Ground parameters at the current frequency (COMMON /GND/).
#[derive(Debug, Clone)]
pub(crate) struct Gnd {
    /// 1 in free space, 2 when a ground plane (and so an image) is present.
    pub ksymp: usize,
    /// Far-field mode of the RP card.
    pub ifar: i32,
    /// 0 reflection coefficients, 1 perfect ground, 2 Sommerfeld-Norton.
    pub iperf: i32,
    /// Number of radials in the ground-screen approximation.
    pub nradl: i32,
    pub t2: f64,
    /// Cliff edge distance and height, in wavelengths.
    pub cl: f64,
    pub ch: f64,
    /// Radial length and radius, in wavelengths.
    pub scrwl: f64,
    pub scrwr: f64,
    /// 1 / sqrt(complex relative permittivity) of the ground.
    pub zrati: Complex64,
    /// The same for the second medium.
    pub zrati2: Complex64,
    pub t1: Complex64,
    /// (eps - 1) / (eps + 1) for the Sommerfeld ground, 1 otherwise.
    pub frati: Complex64,
}

impl Default for Gnd {
    fn default() -> Self {
        Gnd {
            ksymp: 1,
            ifar: -1,
            iperf: 0,
            nradl: 0,
            t2: 0.0,
            cl: 0.0,
            ch: 0.0,
            scrwl: 0.0,
            scrwr: 0.0,
            zrati: C1,
            zrati2: C1,
            t1: C0,
            frati: C1,
        }
    }
}

/// Geometry of one ground-wave evaluation (COMMON /GWAV/).
#[derive(Debug, Clone, Default)]
pub(crate) struct Gwav {
    pub r1: f64,
    pub r2: f64,
    pub zmh: f64,
    pub zph: f64,
    pub u: Complex64,
    pub u2: Complex64,
    pub xx1: Complex64,
    pub xx2: Complex64,
}

/// Observation point of the Sommerfeld field routines (COMMON /INCOM/).
#[derive(Debug, Clone, Default)]
pub(crate) struct Incom {
    pub isnor: i32,
    pub xo: f64,
    pub yo: f64,
    pub zo: f64,
    pub sn: f64,
    pub xsn: f64,
    pub ysn: f64,
}

/// Components of basis functions on segments (COMMON /SEGJ/): segment numbers
/// (1-based) and the constant, sine and cosine coefficients.
#[derive(Debug, Clone, Default)]
pub(crate) struct Basis {
    pub jco: Vec<i32>,
    pub ax: Vec<f64>,
    pub bx: Vec<f64>,
    pub cx: Vec<f64>,
}

impl Basis {
    fn clear(&mut self) {
        self.jco.clear();
        self.ax.clear();
        self.bx.clear();
        self.cx.clear();
    }

    fn push(&mut self, jco: i32, ax: f64, bx: f64, cx: f64) {
        self.jco.push(jco);
        self.ax.push(ax);
        self.bx.push(bx);
        self.cx.push(cx);
    }
}

/// NEC-2 solver state for one structure.
#[derive(Debug, Clone)]
pub struct Engine {
    // --- structure (COMMON /DATA/) ---
    n: usize,
    np: usize,
    /// Number of patches, and of patches in one symmetric cell.
    m: usize,
    mp: usize,
    /// Patches in metres, and in wavelengths at the current frequency.
    patches_m: Vec<Patch>,
    patches: Vec<Patch>,
    ipsym: i32,
    itag: Vec<i32>,
    icon1: Vec<i32>,
    icon2: Vec<i32>,
    /// Segment centres, lengths and radii in metres.
    xm: Vec<f64>,
    ym: Vec<f64>,
    zm: Vec<f64>,
    sim: Vec<f64>,
    bim: Vec<f64>,
    /// The same in wavelengths at the current frequency.
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    si: Vec<f64>,
    bi: Vec<f64>,
    cab: Vec<f64>,
    sab: Vec<f64>,
    salp: Vec<f64>,
    wlam: f64,
    fmhz: f64,

    dj: SourceSeg,
    gnd: Gnd,
    gwav: Gwav,
    incom: Incom,
    segj: Basis,
    sommerfeld: somnec::Tables,

    // --- current expansion (COMMON /CRNT/) ---
    air: Vec<f64>,
    aii: Vec<f64>,
    bir: Vec<f64>,
    bii: Vec<f64>,
    cir: Vec<f64>,
    cii: Vec<f64>,
    cur: Vec<Complex64>,

    // --- interaction matrix ---
    cm: Vec<Complex64>,
    ip: Vec<i32>,
    /// Symmetry transform matrix and the number of symmetric sections.
    ssx: Vec<Complex64>,
    nop: usize,

    // --- sources (COMMON /VSORC/), segment numbers are 1-based ---
    isant: Vec<usize>,
    vsant: Vec<Complex64>,
    ivqd: Vec<usize>,
    vqd: Vec<Complex64>,
    iqds: Vec<usize>,
    vqds: Vec<Complex64>,

    // --- loading (COMMON /ZLOAD/) ---
    nload: usize,
    zarray: Vec<Complex64>,

    // --- networks: solution kept between excitations (COMMON /NETCX/) ---
    net: network::NetState,
    pin: f64,
    pnls: f64,
    zped: Complex64,

    // --- ground parameters as entered (COMMON /SAVE/) ---
    epsr: f64,
    sig: f64,
    scrwlt: f64,
    scrwrt: f64,

    diag: Diag,
    /// Smallest ratio of the smallest to the largest pivot over the blocks
    /// of the last factorization.
    pivot_ratio: f64,
}

impl Engine {
    /// Creates the solver state for a finished structure.
    pub fn new(structure: &Structure) -> Self {
        let n = structure.n();
        Engine {
            n,
            np: structure.np,
            m: structure.m(),
            mp: structure.mp,
            patches_m: structure.patches.clone(),
            patches: structure.patches.clone(),
            ipsym: structure.ipsym,
            itag: structure.itag.clone(),
            icon1: structure.icon1.clone(),
            icon2: structure.icon2.clone(),
            xm: structure.x.clone(),
            ym: structure.y.clone(),
            zm: structure.z.clone(),
            sim: structure.si.clone(),
            bim: structure.bi.clone(),
            x: structure.x.clone(),
            y: structure.y.clone(),
            z: structure.z.clone(),
            si: structure.si.clone(),
            bi: structure.bi.clone(),
            cab: structure.cab.clone(),
            sab: structure.sab.clone(),
            salp: structure.salp.clone(),
            wlam: 1.0,
            fmhz: crate::consts::CVEL,
            dj: SourceSeg { rkh: 1.0, ..SourceSeg::default() },
            gnd: Gnd::default(),
            gwav: Gwav::default(),
            incom: Incom::default(),
            segj: Basis::default(),
            sommerfeld: somnec::Tables::default(),
            air: vec![0.0; n],
            aii: vec![0.0; n],
            bir: vec![0.0; n],
            bii: vec![0.0; n],
            cir: vec![0.0; n],
            cii: vec![0.0; n],
            cur: vec![C0; n + 3 * structure.m()],
            cm: Vec::new(),
            ip: vec![0; n + 2 * structure.m()],
            ssx: Vec::new(),
            nop: 1,
            isant: Vec::new(),
            vsant: Vec::new(),
            ivqd: Vec::new(),
            vqd: Vec::new(),
            iqds: Vec::new(),
            vqds: Vec::new(),
            nload: 0,
            zarray: Vec::new(),
            net: network::NetState::default(),
            pin: 0.0,
            pnls: 0.0,
            zped: C0,
            epsr: 0.0,
            sig: 0.0,
            scrwlt: 0.0,
            scrwrt: 0.0,
            diag: Diag::default(),
            pivot_ratio: 1.0,
        }
    }

    /// Sets the frequency and scales the structure to wavelengths.
    pub fn set_frequency(&mut self, fmhz: f64) {
        use crate::consts::CVEL;
        self.fmhz = fmhz;
        let fr = fmhz / CVEL;
        self.wlam = CVEL / fmhz;
        for i in 0..self.n {
            self.x[i] = self.xm[i] * fr;
            self.y[i] = self.ym[i] * fr;
            self.z[i] = self.zm[i] * fr;
            self.si[i] = self.sim[i] * fr;
            self.bi[i] = self.bim[i] * fr;
        }
        let fr2 = fr * fr;
        for (patch, metric) in self.patches.iter_mut().zip(&self.patches_m) {
            patch.x = metric.x * fr;
            patch.y = metric.y * fr;
            patch.z = metric.z * fr;
            patch.area = metric.area * fr2;
        }
    }

    /// Number of unknowns: one for every segment and two for every patch.
    pub(crate) fn neq(&self) -> usize {
        self.n + 2 * self.m
    }

    /// Number of unknowns in one symmetric cell.
    pub(crate) fn npeq(&self) -> usize {
        self.np + 2 * self.mp
    }

    /// Number of patches.
    pub fn m(&self) -> usize {
        self.m
    }

    /// The patches, with centres in wavelengths and areas in square
    /// wavelengths at the current frequency.
    pub fn patches(&self) -> &[Patch] {
        &self.patches
    }

    /// Surface current of the last solution on every patch, in amperes per
    /// metre: its x, y and z components.
    pub fn patch_currents(&self) -> &[Complex64] {
        &self.cur[self.n..]
    }

    /// Loads patch `j` (0-based) as the source of the patch field routines.
    pub(crate) fn set_source_patch(&mut self, j: usize) {
        let p = self.patches[j];
        self.dj.s = p.area;
        self.dj.xj = p.x;
        self.dj.yj = p.y;
        self.dj.zj = p.z;
        self.dj.t1xj = p.t1x;
        self.dj.t1yj = p.t1y;
        self.dj.t1zj = p.t1z;
        self.dj.t2xj = p.t2x;
        self.dj.t2yj = p.t2y;
        self.dj.t2zj = p.t2z;
    }

    /// Number of segments.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Frequency in MHz.
    pub fn frequency(&self) -> f64 {
        self.fmhz
    }

    /// Wavelength in metres.
    pub fn wavelength(&self) -> f64 {
        self.wlam
    }

    /// Tag of segment `i` (0-based).
    pub fn tag(&self, i: usize) -> i32 {
        self.itag[i]
    }

    /// Centre of segment `i` (0-based) in wavelengths.
    pub fn segment_center(&self, i: usize) -> [f64; 3] {
        [self.x[i], self.y[i], self.z[i]]
    }

    /// Length of segment `i` (0-based) in wavelengths.
    pub fn segment_length(&self, i: usize) -> f64 {
        self.si[i]
    }

    /// Currents of the last solution at the segment centres, in amperes per
    /// wavelength (multiply by [`Engine::wavelength`] for amperes).
    pub fn currents(&self) -> &[Complex64] {
        &self.cur[..self.n]
    }

    /// Charge-density coefficients of the last solution: the sine terms of the
    /// current expansion `(bir, bii)` of segment `i` (0-based).
    pub fn sine_term(&self, i: usize) -> (f64, f64) {
        (self.bir[i], self.bii[i])
    }

    /// Total power delivered by the sources in the last solution, watts.
    pub fn input_power(&self) -> f64 {
        self.pin
    }

    /// Power lost in networks in the last solution, watts.
    pub fn network_loss(&self) -> f64 {
        self.pnls
    }

    /// Input impedance of the last source handled by [`Engine::solve`].
    pub fn source_impedance(&self) -> Complex64 {
        self.zped
    }

    /// Power dissipated in the structure loads for the last solution, watts.
    pub fn structure_loss(&self) -> f64 {
        let mut ploss = 0.0;
        if self.nload != 0 {
            for i in 0..self.n {
                let curi = self.cur[i] * self.wlam;
                let cmag = curi.norm();
                if self.zarray[i].re.abs() >= 1.0e-20 {
                    ploss += 0.5 * cmag * cmag * self.zarray[i].re * self.si[i];
                }
            }
        }
        ploss
    }

    /// Returns and clears the warnings collected by the numerical routines
    /// (for example a pivot that is nearly zero).
    pub fn take_warnings(&mut self) -> Vec<String> {
        std::mem::take(&mut self.diag.warnings)
    }

    /// Ratio of the smallest to the largest pivot of the last factorization,
    /// the minimum over the blocks of a symmetric structure. It is a cheap
    /// indicator of the conditioning of the interaction matrix: about
    /// `-log10` of it decimal digits of the solution are lost. Sound models
    /// give more than 0.01; below 1e-5 a warning is issued.
    pub fn pivot_ratio(&self) -> f64 {
        self.pivot_ratio
    }

    /// Turns a failure reported by an integration routine into an error that
    /// names the source and the observation segment (numbered from 1).
    pub(crate) fn check_integration(&mut self, source: usize, observation: Option<usize>) -> Result<()> {
        match self.diag.fault.take() {
            None => Ok(()),
            Some(routine) => Err(Error::Integration { routine, source: Some(source), observation }),
        }
    }
}
