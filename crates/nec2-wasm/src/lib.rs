//! NEC-2 for JavaScript: runs a deck of NEC cards with the `nec2` engine and
//! returns the report and the main results as numbers.
//!
//! The results are plain JavaScript objects. A pattern is stored by columns,
//! one array for every quantity, because a pattern of the whole sphere has
//! tens of thousands of directions.

use nec2_deck::{JobResults, Options, Output, OverlapCheck, PortRecord, RadiationPattern, Solution};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[wasm_bindgen(typescript_custom_section)]
const TYPES: &str = r#"
/** A complex number: real and imaginary part. */
export type Complex = [number, number];

export interface RunOptions {
  /** Refuse a structure with more equations (one for a segment, two for a patch). */
  maxEquations?: number;
  /** What to do when wires occupy the same space. The default is "warn". */
  overlap?: "warn" | "error" | "off";
}

export interface Output {
  /** The report, in the layout of the NEC-2 output file. */
  report: string;
  /** The message of the error that stopped the run, or null. */
  error: string | null;
  /** The warnings of the report. */
  warnings: string[];
  /** One entry for every job of the deck: one, unless the deck uses NX cards. */
  jobs: Job[];
}

export interface Job {
  segments: Segment[];
  /** One solution for every frequency of a sweep and every direction of an incident wave. */
  solutions: Solution[];
}

/** A segment of the structure: coordinates and radius in metres. */
export interface Segment {
  tag: number;
  end1: [number, number, number];
  end2: [number, number, number];
  radius: number;
}

export interface Solution {
  /** Frequency in MHz. */
  frequency: number;
  /** For a plane wave: theta, phi and the polarization angle eta, in degrees. */
  incidentWave: [number, number, number] | null;
  sources: Port[];
  networkPorts: Port[];
  /** Current at the centre of every segment, in amperes. */
  currents: Complex[];
  /** Power in watts, for voltage sources. */
  power: { input: number; radiated: number; structureLoss: number; networkLoss: number } | null;
  patterns: Pattern[];
}

/** Voltage, current and power at a source or at a port of a network. */
export interface Port {
  tag: number;
  /** Number of the segment, from 1. */
  segment: number;
  voltage: Complex;
  current: Complex;
  impedance: Complex;
  admittance: Complex;
  /** Watts. */
  power: number;
}

/**
 * A far-field pattern by columns: entry i of every array belongs to
 * direction i. Angles are in degrees, gains in dB (-999.99 where there is no
 * field), fields in volts per metre.
 */
export interface Pattern {
  theta: number[];
  phi: number[];
  gainTotal: number[];
  gainVertical: number[];
  gainHorizontal: number[];
  gainMajor: number[];
  gainMinor: number[];
  axialRatio: number[];
  tilt: number[];
  /** "", "LINEAR", "RIGHT" or "LEFT". */
  sense: string[];
  eThetaMagnitude: number[];
  eThetaPhase: number[];
  ePhiMagnitude: number[];
  ePhiPhase: number[];
  /** The average power gain, if the RP card asked for it and it is defined. */
  averageGain: number | null;
}
"#;

#[derive(Deserialize, Default)]
#[serde(rename_all = "camelCase")]
struct RunOptions {
    max_equations: Option<usize>,
    overlap: Option<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct OutputView<'a> {
    report: &'a str,
    error: Option<&'a str>,
    warnings: &'a [String],
    jobs: Vec<JobView>,
}

#[derive(Serialize)]
struct JobView {
    segments: Vec<SegmentView>,
    solutions: Vec<SolutionView>,
}

#[derive(Serialize)]
struct SegmentView {
    tag: i32,
    end1: [f64; 3],
    end2: [f64; 3],
    radius: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct SolutionView {
    frequency: f64,
    incident_wave: Option<[f64; 3]>,
    sources: Vec<PortView>,
    network_ports: Vec<PortView>,
    currents: Vec<[f64; 2]>,
    power: Option<PowerView>,
    patterns: Vec<PatternView>,
}

#[derive(Serialize)]
struct PortView {
    tag: i32,
    segment: usize,
    voltage: [f64; 2],
    current: [f64; 2],
    impedance: [f64; 2],
    admittance: [f64; 2],
    power: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PowerView {
    input: f64,
    radiated: f64,
    structure_loss: f64,
    network_loss: f64,
}

#[derive(Serialize, Default)]
#[serde(rename_all = "camelCase")]
struct PatternView {
    theta: Vec<f64>,
    phi: Vec<f64>,
    gain_total: Vec<f64>,
    gain_vertical: Vec<f64>,
    gain_horizontal: Vec<f64>,
    gain_major: Vec<f64>,
    gain_minor: Vec<f64>,
    axial_ratio: Vec<f64>,
    tilt: Vec<f64>,
    sense: Vec<&'static str>,
    e_theta_magnitude: Vec<f64>,
    e_theta_phase: Vec<f64>,
    e_phi_magnitude: Vec<f64>,
    e_phi_phase: Vec<f64>,
    average_gain: Option<f64>,
}

fn port(record: &PortRecord) -> PortView {
    PortView {
        tag: record.tag,
        segment: record.segment,
        voltage: [record.voltage.re, record.voltage.im],
        current: [record.current.re, record.current.im],
        impedance: [record.impedance.re, record.impedance.im],
        admittance: [record.admittance.re, record.admittance.im],
        power: record.power,
    }
}

fn pattern(pattern: &RadiationPattern) -> PatternView {
    let mut view =
        PatternView { average_gain: pattern.average_gain.and_then(|average| average.gain), ..Default::default() };
    for row in &pattern.rows {
        view.theta.push(row.theta);
        view.phi.push(row.phi);
        view.gain_total.push(row.gain_total);
        view.gain_vertical.push(row.gain_vertical);
        view.gain_horizontal.push(row.gain_horizontal);
        view.gain_major.push(row.gain_major);
        view.gain_minor.push(row.gain_minor);
        view.axial_ratio.push(row.axial_ratio);
        view.tilt.push(row.tilt);
        view.sense.push(row.sense.trim());
        view.e_theta_magnitude.push(row.e_theta_mag);
        view.e_theta_phase.push(row.e_theta_phase);
        view.e_phi_magnitude.push(row.e_phi_mag);
        view.e_phi_phase.push(row.e_phi_phase);
    }
    view
}

fn solution(solution: &Solution) -> SolutionView {
    SolutionView {
        frequency: solution.frequency,
        incident_wave: solution.incident_wave,
        sources: solution.sources.iter().map(port).collect(),
        network_ports: solution.network_ports.iter().map(port).collect(),
        currents: solution.currents.iter().map(|current| [current.re, current.im]).collect(),
        power: solution.power.map(|power| PowerView {
            input: power.input,
            radiated: power.radiated,
            structure_loss: power.structure_loss,
            network_loss: power.network_loss,
        }),
        patterns: solution.patterns.iter().map(pattern).collect(),
    }
}

fn job(job: &JobResults) -> JobView {
    let segment = |segment: &nec2_deck::Segment| SegmentView {
        tag: segment.tag,
        end1: segment.end1,
        end2: segment.end2,
        radius: segment.radius,
    };
    JobView {
        segments: job.segments.iter().map(segment).collect(),
        solutions: job.solutions.iter().map(solution).collect(),
    }
}

fn view(output: &Output) -> OutputView<'_> {
    OutputView {
        report: &output.report,
        error: output.error.as_deref(),
        warnings: &output.warnings,
        jobs: output.jobs.iter().map(job).collect(),
    }
}

fn to_js(output: &Output) -> Result<JsValue, JsError> {
    // Plain objects and `null` for what is absent, as JSON would have them.
    let serializer = serde_wasm_bindgen::Serializer::json_compatible();
    Ok(view(output).serialize(&serializer)?)
}

/// Runs a deck of NEC cards. An error of the deck does not throw: it is in
/// `error` of the result, with the report up to it.
#[wasm_bindgen(unchecked_return_type = "Output")]
pub fn run(deck: &str) -> Result<JsValue, JsError> {
    to_js(&nec2_deck::run(deck))
}

/// Runs a deck of NEC cards with options.
#[wasm_bindgen(js_name = runWith, unchecked_return_type = "Output")]
pub fn run_with(
    deck: &str,
    #[wasm_bindgen(unchecked_param_type = "RunOptions")] options: JsValue,
) -> Result<JsValue, JsError> {
    let given: RunOptions = serde_wasm_bindgen::from_value(options)?;
    let overlap = match given.overlap.as_deref() {
        None | Some("warn") => OverlapCheck::Warn,
        Some("error") => OverlapCheck::Error,
        Some("off") => OverlapCheck::Off,
        Some(other) => return Err(JsError::new(&format!("overlap must be warn, error or off, not {other}"))),
    };
    let options = Options { overlap, max_equations: given.max_equations };
    to_js(&nec2_deck::run_with(deck, &options))
}

/// The version of the package.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_view_holds_the_numbers_of_the_run() {
        let deck = "CM Dipole\nCE\nGW 1 9 0 0 -0.25 0 0 0.25 0.001\nGE 0\nEX 0 1 5 0 1 0\nFR 0 2 0 0 290 10\n\
                    RP 0 3 2 1000 0 0 45 90\nEN\n";
        let output = nec2_deck::run(deck);
        let view = view(&output);
        assert_eq!(view.error, None);
        let [job] = &view.jobs[..] else { panic!("one job") };
        assert_eq!((job.segments.len(), job.solutions.len()), (9, 2));
        let solution = &job.solutions[1];
        assert_eq!(solution.frequency, 300.0);
        assert_eq!((solution.currents.len(), solution.sources.len()), (9, 1));
        let impedance = output.jobs[0].solutions[1].sources[0].impedance;
        assert_eq!(solution.sources[0].impedance, [impedance.re, impedance.im]);
        let [pattern] = &solution.patterns[..] else { panic!("one pattern") };
        assert_eq!(pattern.theta, [0.0, 45.0, 90.0, 0.0, 45.0, 90.0]);
        assert_eq!(pattern.sense[1], "LINEAR");
        assert!(pattern.gain_total[2] > 2.0);
    }
}
