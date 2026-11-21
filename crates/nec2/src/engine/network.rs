//! Solution for the currents, including non-radiating networks and
//! transmission lines (NT and TL cards).
//!
//! Ported from `netwk` (`network.c`) of nec2c 1.3.1. Three things follow the
//! original NEC-2 Fortran instead, because nec2c 1.3.1 departs from it (see
//! `docs/nec2c-deviations.md`): the network solution is kept between
//! excitations, a source that is not on a network port is not counted as
//! network loss, and the admittance-asymmetry check works without networks.

use num_complex::Complex64;

use super::matrix::{factr, solve};
use super::{C0, C1, Engine, all_finite};
use crate::consts::{TP, cmplx};
use crate::error::{Error, Result};

/// Kind of a two-port network.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NetworkKind {
    /// General two-port given by its short-circuit admittance matrix (NT).
    TwoPort,
    /// Transmission line (TL).
    Line,
    /// Transmission line with a 180 degree phase reversal (TL, negative impedance).
    CrossedLine,
}

/// A two-port network between two segments.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Network {
    pub kind: NetworkKind,
    /// Segments (1-based) that ports 1 and 2 are connected to.
    pub iseg1: usize,
    pub iseg2: usize,
    /// For a two-port: Y11, Y12 and Y22 in siemens (real, imaginary).
    /// For a line: `x11r` is the characteristic impedance in ohms, `x11i`
    /// the length in metres, `x12r`/`x12i` the shunt admittance at end 1
    /// and `x22r`/`x22i` the shunt admittance at end 2.
    pub x11r: f64,
    pub x11i: f64,
    pub x12r: f64,
    pub x12i: f64,
    pub x22r: f64,
    pub x22i: f64,
}

/// Voltage, current and power at a source or at a network port.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PortRecord {
    pub tag: i32,
    /// Segment number, 1-based.
    pub segment: usize,
    pub voltage: Complex64,
    pub current: Complex64,
    pub impedance: Complex64,
    pub admittance: Complex64,
    pub power: f64,
}

/// Relative asymmetry of the driving-point admittance matrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Asymmetry {
    pub max: f64,
    /// The pair of segments (1-based) with the largest asymmetry.
    pub segment1: usize,
    pub segment2: usize,
    pub rms: f64,
}

/// Results of [`Engine::solve`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SolveReport {
    pub asymmetry: Option<Asymmetry>,
    /// Structure excitation at the network connection points.
    pub network_ports: Vec<PortRecord>,
    /// Antenna input parameters, one record per source.
    pub sources: Vec<PortRecord>,
}

/// Network equations of the last solution, kept so that further excitations
/// can reuse them.
#[derive(Debug, Clone, Default)]
pub(crate) struct NetState {
    nteq: usize,
    ntsc: usize,
    ndimn: usize,
    cmn: Vec<Complex64>,
    ipnt: Vec<i32>,
    /// Segments (1-based) of the network ports without a source.
    nteqa: Vec<usize>,
    /// Segments (1-based) of the network ports that carry a voltage source.
    ntsca: Vec<usize>,
    vsrc: Vec<Complex64>,
    rhnx: Vec<Complex64>,
}

impl Engine {
    fn port_record(&self, segment: usize, voltage: Complex64, current: Complex64) -> PortRecord {
        PortRecord {
            tag: self.itag[segment],
            segment: segment + 1,
            voltage,
            current,
            impedance: voltage / current,
            admittance: current / voltage,
            power: 0.5 * (voltage * current.conj()).re,
        }
    }

    /// Solves for the structure currents for the excitation set by
    /// [`Engine::excite`], including the effect of `networks`.
    ///
    /// With `asymmetry` the relative asymmetry of the driving-point admittance
    /// matrix is computed as a check of the solution. `reuse` keeps the
    /// network equations of the previous call, which is valid as long as the
    /// matrix, the networks and the voltage sources are unchanged (for
    /// example between the incidence angles of a plane-wave sweep).
    pub fn solve(&mut self, networks: &[Network], asymmetry: bool, reuse: bool) -> Result<SolveReport> {
        let wlam = self.wlam;
        let nonet = networks.len();
        let nsant = self.isant.len();
        let mut report = SolveReport::default();
        let mut einc = std::mem::take(&mut self.cur);
        // Room for the three current components that `cabc` writes for
        // every patch.
        let mut rhs = vec![C0; einc.len()];
        self.pin = 0.0;
        self.pnls = 0.0;

        if !reuse {
            let ndimn = 2 * nonet + nsant;
            let mut st = NetState {
                ndimn,
                cmn: vec![C0; ndimn * ndimn],
                ipnt: vec![0; ndimn],
                rhnx: vec![C0; ndimn],
                ..NetState::default()
            };

            // Relative asymmetry of the admittance matrix seen from the
            // network ports and sources.
            if asymmetry {
                let mut ports: Vec<usize> = Vec::new();
                let candidates =
                    networks.iter().flat_map(|net| [net.iseg1, net.iseg2]).chain(self.isant.iter().copied());
                for segment in candidates {
                    if !ports.contains(&segment) {
                        ports.push(segment);
                    }
                }
                let irow1 = ports.len();
                if irow1 >= 2 {
                    let mut ymat = vec![C0; irow1 * irow1];
                    for i in 0..irow1 {
                        let isc1 = ports[i] - 1;
                        let asmx = self.si[isc1];
                        rhs.fill(C0);
                        rhs[isc1] = C1;
                        self.solves(&mut rhs);
                        self.cabc(&mut rhs)?;
                        for j in 0..irow1 {
                            ymat[j + i * irow1] = rhs[ports[j] - 1] / asmx;
                        }
                    }

                    let mut asmx = 0.0;
                    let mut asa = 0.0;
                    let (mut segment1, mut segment2) = (0, 0);
                    for i in 1..irow1 {
                        for j in 0..i {
                            let cux = ymat[i + j * irow1];
                            let pwr = ((cux - ymat[j + i * irow1]) / cux).norm();
                            asa += pwr * pwr;
                            if pwr < asmx {
                                continue;
                            }
                            asmx = pwr;
                            segment1 = ports[i];
                            segment2 = ports[j];
                        }
                    }
                    let rms = (asa * 2.0 / (irow1 * (irow1 - 1)) as f64).sqrt();
                    report.asymmetry = Some(Asymmetry { max: asmx, segment1, segment2, rms });
                }
            }

            // Network equations: sort the ports into those with and without
            // a voltage source and fill the short-circuit admittances.
            if nonet != 0 {
                for net in networks {
                    let (nseg1, nseg2) = (net.iseg1, net.iseg2);
                    let (y11, y12, y22);
                    if net.kind == NetworkKind::TwoPort {
                        y11 = cmplx(net.x11r, net.x11i);
                        y12 = cmplx(net.x12r, net.x12i);
                        y22 = cmplx(net.x22r, net.x22i);
                    } else {
                        // Lossless line of electrical length `bl` with shunt
                        // admittances at its ends.
                        let bl = TP * net.x11i / wlam;
                        let mut y12i = 1.0 / (net.x11r * bl.sin());
                        let y11i = -y12i * bl.cos();
                        y11 = cmplx(net.x12r, y11i + net.x12i);
                        y22 = cmplx(net.x22r, y11i + net.x22i);
                        if net.kind == NetworkKind::CrossedLine {
                            y12i = -y12i;
                        }
                        y12 = cmplx(0.0, y12i);
                    }

                    let isc1 = self.isant.iter().position(|&s| s == nseg1);
                    let irow1 = st.assign_row(nseg1, isc1.map(|i| self.vsant[i]));
                    let isc2 = self.isant.iter().position(|&s| s == nseg2);
                    let irow2 = st.assign_row(nseg2, isc2.map(|i| self.vsant[i]));
                    let ndimn = st.ndimn;

                    match isc1 {
                        None => {
                            st.cmn[irow1 + irow1 * ndimn] -= y11 * self.si[nseg1 - 1];
                            st.cmn[irow1 + irow2 * ndimn] -= y12 * self.si[nseg1 - 1];
                        }
                        Some(isc) => {
                            st.rhnx[irow1] += y11 * self.vsant[isc] / wlam;
                            st.rhnx[irow2] += y12 * self.vsant[isc] / wlam;
                        }
                    }
                    match isc2 {
                        None => {
                            st.cmn[irow2 + irow2 * ndimn] -= y22 * self.si[nseg2 - 1];
                            st.cmn[irow2 + irow1 * ndimn] -= y12 * self.si[nseg2 - 1];
                        }
                        Some(isc) => {
                            st.rhnx[irow1] += y12 * self.vsant[isc] / wlam;
                            st.rhnx[irow2] += y22 * self.vsant[isc] / wlam;
                        }
                    }
                }

                // Add the admittances of the structure between the ports.
                for i in 0..st.nteq {
                    rhs.fill(C0);
                    rhs[st.nteqa[i] - 1] = C1;
                    self.solves(&mut rhs);
                    self.cabc(&mut rhs)?;
                    for j in 0..st.nteq {
                        st.cmn[i + j * st.ndimn] += rhs[st.nteqa[j] - 1];
                    }
                }

                if let Err(equation) = factr(st.nteq, &mut st.cmn, &mut st.ipnt, st.ndimn, &mut self.diag.warnings) {
                    self.cur = einc;
                    return Err(Error::SingularMatrix { matrix: "NETWORK", equation, frequency: self.fmhz });
                }
            }
            self.net = st;
        }

        let mut rhnt = vec![C0; self.net.ndimn];
        let ntsc = if nonet != 0 {
            // Add to the network right-hand side the currents that the
            // excitation alone drives through the ports.
            rhs.copy_from_slice(&einc);
            self.solves(&mut rhs);
            self.cabc(&mut rhs)?;
            let st = &self.net;
            for i in 0..st.nteq {
                rhnt[i] = st.rhnx[i] + rhs[st.nteqa[i] - 1];
            }

            // Solve the network equations for the port voltages, apply them
            // to the structure and solve for the currents.
            solve(st.nteq, &st.cmn, &st.ipnt, &mut rhnt, st.ndimn);
            for i in 0..st.nteq {
                einc[st.nteqa[i] - 1] -= rhnt[i];
            }
            self.solves(&mut einc);
            self.cabc(&mut einc)?;

            let st = &self.net;
            for i in 0..st.nteq {
                let irow1 = st.nteqa[i] - 1;
                let vlt = rhnt[i] * self.si[irow1] * wlam;
                let record = self.port_record(irow1, vlt, einc[irow1] * wlam);
                self.zped = record.impedance;
                self.pnls -= record.power;
                report.network_ports.push(record);
            }
            for i in 0..st.ntsc {
                let irow1 = st.ntsca[i] - 1;
                let record = self.port_record(irow1, st.vsrc[i], einc[irow1] * wlam);
                self.zped = record.impedance;
                self.pnls -= record.power;
                report.network_ports.push(record);
            }
            st.ntsc
        } else {
            // No networks.
            self.solves(&mut einc);
            self.cabc(&mut einc)?;
            0
        };

        // Antenna input parameters.
        for i in 0..nsant {
            let isc1 = self.isant[i] - 1;
            let vlt = self.vsant[i];
            let st = &self.net;
            let port = if ntsc == 0 { None } else { st.ntsca.iter().position(|&s| s == isc1 + 1) };
            let cux = match port {
                None => einc[isc1] * wlam,
                Some(j) => {
                    let irow1 = st.ndimn - (j + 1);
                    let mut cux = st.rhnx[irow1];
                    for k in 0..st.nteq {
                        cux -= st.cmn[k + irow1 * st.ndimn] * rhnt[k];
                    }
                    (einc[isc1] + cux) * wlam
                }
            };
            let record = self.port_record(isc1, vlt, cux);
            self.zped = record.impedance;
            self.pin += record.power;
            if port.is_some() {
                self.pnls += record.power;
            }
            report.sources.push(record);
        }

        for i in 0..self.ivqd.len() {
            let isc1 = self.ivqd[i] - 1;
            let vlt = self.vqd[i];
            let a = cmplx(self.air[isc1], self.aii[isc1]);
            let b = cmplx(self.bir[isc1], self.bii[isc1]);
            let c = cmplx(self.cir[isc1], self.cii[isc1]);
            let half = self.si[isc1] * TP * 0.5;
            let cux = (a - b * half.sin() + c * half.cos()) * wlam;
            let record = self.port_record(isc1, vlt, cux);
            self.zped = record.impedance;
            self.pin += record.power;
            report.sources.push(record);
        }

        // A solution that is not finite must not reach the results.
        let finite = all_finite(&einc);
        self.cur = einc;
        if !finite {
            return Err(Error::NonFiniteSolution { frequency: self.fmhz });
        }
        Ok(report)
    }
}

impl NetState {
    /// Returns the row of the network equations for the port on `segment`.
    /// Ports without a source are numbered from the top, ports with a voltage
    /// source (`source` holds its voltage) from the bottom.
    fn assign_row(&mut self, segment: usize, source: Option<Complex64>) -> usize {
        match source {
            None => match self.nteqa.iter().position(|&s| s == segment) {
                Some(i) => i,
                None => {
                    self.nteqa.push(segment);
                    self.nteq += 1;
                    self.nteq - 1
                }
            },
            Some(voltage) => match self.ntsca.iter().position(|&s| s == segment) {
                Some(i) => self.ndimn - (i + 1),
                None => {
                    self.ntsca.push(segment);
                    self.vsrc.push(voltage);
                    self.ntsc += 1;
                    self.ndimn - self.ntsc
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::Excitation;
    use crate::geometry::Geometry;

    #[test]
    fn a_solution_that_is_not_finite_is_an_error() {
        let mut geometry = Geometry::new();
        geometry.wire(1, 9, [0.0, 0.0, -0.25], [0.0, 0.0, 0.25], 0.001, 1.0, 1.0).unwrap();
        let mut engine = Engine::new(&geometry.connect(0).unwrap());
        engine.set_frequency(299.8);
        engine.fill_and_factor(1.0, false).unwrap();
        assert!(engine.pivot_ratio() > 0.05 && engine.take_warnings().is_empty());

        engine.set_sources(&[(5, cmplx(f64::NAN, 0.0))], &[]);
        engine.excite(&Excitation::Sources).unwrap();
        assert_eq!(engine.solve(&[], false, false).err(), Some(Error::NonFiniteSolution { frequency: 299.8 }));
        // The engine stays usable.
        assert_eq!(engine.currents().len(), 9);
        engine.set_sources(&[(5, C1)], &[]);
        engine.excite(&Excitation::Sources).unwrap();
        assert!(engine.solve(&[], false, false).is_ok());
    }
}
