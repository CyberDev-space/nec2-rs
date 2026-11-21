//! Runs the decks in `tests/fixtures/decks` and compares the reports with
//! the reference reports of nec2c 1.3.1 in `tests/fixtures/nec2c`.

mod support;

use std::path::PathBuf;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures")
}

/// Decks whose results nec2c 1.3.1 gets slightly wrong; they are compared
/// with the reports of nec2c with the fixes of `scripts/patch-nec2c.py`
/// (see `docs/nec2c-deviations.md`).
const FIXED_REFERENCE: &[&str] = &[
    "s4_som_gp_radials",
    "s4_som_hdipole_high",
    "s4_som_hdipole_low",
    "s4_som_hf_dipole",
    "s4_som_sweep",
    "s4_som_vertical",
    "s6_asymmetry",
    "s6_ld_wire_hf",
    "s6_ld_wire_thin",
    "s6_rp_range",
    "s6_som_long_wire",
    "s6_som_sea_array",
    "s6_tl_two_sources",
    "s8_patch_current_on_patch",
    "s8_patch_strip",
    "s8_patch_stub_first",
];

fn read(directory: &str, name: &str, extension: &str) -> String {
    let path = fixtures().join(format!("{directory}/{name}.{extension}"));
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn load(name: &str) -> (String, String) {
    let directory = if FIXED_REFERENCE.contains(&name) { "nec2c-fixed" } else { "nec2c" };
    (read("decks", name, "nec"), read(directory, name, "out"))
}

/// Runs a deck and compares the report, without the warnings that nec2-rs
/// adds, with the reference. Returns the report.
fn check(name: &str, units: Option<f64>) -> String {
    let (deck, reference) = load(name);
    let output = nec2_deck::run(&deck);
    assert_eq!(output.error, None, "{name}: the run stopped with an error");
    let mismatches = support::compare(&reference, &support::without_warnings(&output.report), units, false);
    assert!(mismatches.is_empty(), "{}", support::describe(name, &mismatches));
    output.report
}

macro_rules! decks {
    ($($name:ident),* $(,)?) => {
        $(
            #[test]
            fn $name() {
                check(stringify!($name), None);
            }
        )*
    };
}

// Wires in free space.
decks!(
    s1_dipole,
    s1_dipole_ek,
    s1_dipole_msweep,
    s1_dipole_sweep,
    s1_ground_plane,
    s1_inverted_v,
    s1_loop,
    s1_stepped_radius,
    s1_taper,
    s1_two_sources,
    s1_yagi3,
);

// Loads, transmission lines and networks.
decks!(s2_ld_cancel, s2_ld_distributed, s2_ld_lumped, s2_nt_network, s2_nt_tl_mixed, s2_tl_parasitic, s2_tl_phased,);

// Perfect ground and reflection-coefficient ground.
decks!(
    s3_pg_hdipole,
    s3_pg_low_dipole,
    s3_pg_monopole,
    s3_pg_sweep,
    s3_rc_cliff,
    s3_rc_gd,
    s3_rc_ground_wave,
    s3_rc_hdipole,
    s3_rc_monopole,
    s3_rc_monopole_open,
    s3_rc_radials,
    s3_rc_vdipole,
);

// Sommerfeld-Norton ground.
decks!(
    s4_som_gp_radials,
    s4_som_hdipole_high,
    s4_som_hdipole_low,
    s4_som_hf_dipole,
    s4_som_monopole,
    s4_som_sweep,
    s4_som_vertical,
);

// Other geometry cards, sources and print options.
decks!(
    s5_charge_print,
    s5_coupling,
    s5_ex_current_slope,
    s5_ex_elementary,
    s5_ex_elliptic,
    s5_ex_plane_wave,
    s5_geo_arc,
    s5_geo_gm,
    s5_geo_gr,
    s5_geo_gs,
    s5_geo_gx,
    s5_geo_gx3,
    s5_geo_helix,
    s5_kh,
    s5_near_fields,
    s5_normalized_gain,
    s5_two_jobs,
);

// Edge cases, and the places where nec2-rs follows the NEC-2 Fortran
// instead of nec2c.
decks!(
    s6_asymmetry,
    s6_ek_structures,
    s6_elliptic_helix,
    s6_ex5_ek,
    s6_gm_partial,
    s6_gn_cancel,
    s6_gr_loads,
    s6_gx_finite_ground,
    s6_gx_ground,
    s6_iped_max,
    s6_ld_wire_hf,
    s6_ld_wire_thin,
    s6_near_ek,
    s6_plane_wave_ground,
    s6_pt_modes,
    s6_rp_range,
    s6_som_long_wire,
    s6_som_sea_array,
    s6_spiral,
    s6_tl_cancel,
    s6_tl_two_sources,
);

// Surface patches.
decks!(
    s8_patch_bridge,
    s8_patch_coupling,
    s8_patch_current_element,
    s8_patch_current_on_patch,
    s8_patch_cylinder,
    s8_patch_dipole_plate,
    s8_patch_disc_gr,
    s8_patch_ek_slope,
    s8_patch_finite_ground,
    s8_patch_ge1,
    s8_patch_gm_gs,
    s8_patch_gr_wires,
    s8_patch_ground,
    s8_patch_gx_box,
    s8_patch_gx_ground,
    s8_patch_monopole_plate,
    s8_patch_network,
    s8_patch_only_fields,
    s8_patch_only_sweep,
    s8_patch_plate_wave,
    s8_patch_print,
    s8_patch_shapes,
    s8_patch_sommerfeld,
    s8_patch_strip,
    s8_patch_stub_first,
    s8_patch_stub_last,
    s8_patch_sweep_load,
    s8_patch_tri_quad_wave,
    s8_patch_wave_ground,
    s8_patch_wire_end2,
);

/// Decks with PL cards: the plot data is compared with the plot file of
/// nec2c as well.
fn check_plot(name: &str) {
    let (deck, _) = load(name);
    check(name, None);
    let plot = nec2_deck::run(&deck).plot.unwrap_or_else(|| panic!("{name}: no plot data"));
    let mismatches = support::compare(&read("nec2c", name, "plt"), &plot, None, false);
    assert!(mismatches.is_empty(), "{}", support::describe(&format!("{name} (plot file)"), &mismatches));
}

macro_rules! plots {
    ($($name:ident),* $(,)?) => {
        $(
            #[test]
            fn $name() {
                check_plot(stringify!($name));
            }
        )*
    };
}

plots!(s9_plot_currents, s9_plot_near_fields, s9_plot_patches, s9_plot_pattern);

/// A deck without a PL card has no plot data.
#[test]
fn no_plot_data_without_a_pl_card() {
    assert_eq!(nec2_deck::run(&load("s1_dipole").0).plot, None);
}

/// Decks that NEC-2 solves although a card or the model is out of the
/// ordinary. nec2-rs gives the same numbers and adds a warning where the
/// result deserves no trust; `warning` is a part of its text, or empty if
/// there is nothing to warn about.
fn compatible(name: &str, warning: &str) {
    let report = check(name, None);
    let warnings = support::warnings(&report);
    if warning.is_empty() {
        assert!(warnings.is_empty(), "{name}: unexpected {warnings:?}");
    } else {
        assert!(warnings.iter().any(|line| line.contains(warning)), "{name}: {warnings:?}");
    }
}

macro_rules! compatible {
    ($($name:ident => $warning:expr),* $(,)?) => {
        $(
            #[test]
            fn $name() {
                compatible(stringify!($name), $warning);
            }
        )*
    };
}

compatible!(
    // A ground flag and a ground card that do not fit each other.
    s7_ge1_free_space => "",
    s7_ge1_free_space_touch => "SEGMENT ENDS JOINED TO THEIR IMAGE BY THE GE CARD: 1,",
    s7_below_ground_ge0 => "REACH BELOW THE GROUND PLANE OR LIE IN IT: 2, THE FIRST IS SEGMENT 1",
    s7_below_ground_gem1 => "REACH BELOW THE GROUND PLANE OR LIE IN IT: 2,",
    s7_in_ground_plane => "REACH BELOW THE GROUND PLANE OR LIE IN IT: 9,",
    s7_touch_ground_ge0 => "",
    // Counts and tags.
    s7_zero_segments => "HAS NO SEGMENTS AND IS IGNORED",
    s7_negative_tag => "",
    s7_gm_negative_increment => "",
    s7_gx_negative_increment => "",
    s7_gr_negative_increment => "",
    // Ranges that reach beyond the structure.
    s7_ld_past_end => "",
    s7_ld_none => "SEGMENTS 12 TO 14 ARE NOT AMONG THE 9 SEGMENTS WITH TAG 1",
    s7_ld_abs_past_end => "",
    s7_pt_past_end => "",
    s7_pt_none => "SEGMENTS 50 TO 60 ARE NOT AMONG THE 9 SEGMENTS",
    s7_pt_receive_past_end => "",
    // Ground types and pattern modes outside their definition.
    s7_gn_type_3 => "GROUND TYPE 3 IS NOT DEFINED",
    s7_gn_type_minus_2 => "GROUND TYPE -2 IS NOT DEFINED",
    s7_rp4_no_screen => "COMPUTED WITHOUT A SCREEN",
    s7_rp4_no_screen_perfect => "COMPUTED WITHOUT A SCREEN",
    s7_rp56_no_screen => "COMPUTED WITHOUT A SCREEN",
    // Wires in the same space.
    s7_overlap_crossing => "INSIDE THE VOLUME",
    s7_overlap_acute => "SEGMENT 12 (TAG 2) HAS ITS CENTRE INSIDE THE VOLUME OF SEGMENT 4 (TAG 1)",
    // The left-handed helix of NEC-2.
    s7_left_helix_feed => "LEFT-HANDED HELIX",
    // A plane wave from one direction: later requests use the same currents.
    s7_wave_rp_rp => "",
    s7_wave_xq_xq => "",
    s7_wave_xq_rp => "",
    // Patches.
    s7_patch_no_area => "PATCH 1 HAS NO AREA",
    s7_cp_no_segment => "NO SEGMENT 1 OF TAG 0 FOR THE COUPLING",
    // Values that NEC-2 does not define but handles.
    s7_gx_negative_flags => "THE REFLECTION FLAGS -110 ARE NEGATIVE",
    s7_rp_negative_mode => "PATTERN MODE -2 IS NOT DEFINED - AS IN NEC-2, IT IS TAKEN AS MODE 0",
    s7_fr_stepping_type => "STEPPING TYPE 2 IS NOT DEFINED - AS IN NEC-2, THE STEP IS ADDED",
    // The message about a small pivot of the network matrix has its place in the report.
    s7_nt_small_pivot => "",
);

/// A deck that must stop with an error that contains `message`. The error
/// text is this program's own, so the report is not compared with nec2c.
fn check_error(name: &str, message: &str) {
    let (deck, _) = load(name);
    let error = nec2_deck::run(&deck).error.unwrap_or_else(|| panic!("{name} must stop with an error"));
    assert!(error.contains(message), "{name}: {error}");
}

#[test]
fn err_unknown_card() {
    check_error("err_unknown_card", "ERROR AT LINE 6: ZZ 1 2 3\n  UNKNOWN CARD");
}

#[test]
fn err_bad_tag() {
    check_error("err_bad_tag", "NO SEGMENT HAS AN ITAG OF 7");
}

#[test]
fn err_below_ground() {
    check_error("err_below_ground", "GEOMETRY DATA ERROR -- SEGMENT 1 EXTENDS BELOW GROUND");
}

/// Outside the decks listed in `FIXED_REFERENCE` the fixes do not change
/// what nec2c prints, so those decks are compared with the stock program.
#[test]
fn fixes_change_only_the_listed_decks() {
    let ignore_run_times =
        |text: &str| text.lines().filter(|line| !line.contains("msec")).collect::<Vec<_>>().join("\n");
    for entry in std::fs::read_dir(fixtures().join("decks")).expect("deck directory") {
        let file_name = entry.expect("directory entry").file_name();
        let name = file_name.to_str().and_then(|name| name.strip_suffix(".nec")).expect("deck file name");
        let same = ignore_run_times(&read("nec2c", name, "out")) == ignore_run_times(&read("nec2c-fixed", name, "out"));
        assert_eq!(same, !FIXED_REFERENCE.contains(&name), "{name}");
    }
}

/// Prints how every deck compares, strictly and with the noise rules.
/// Run with `cargo test -p nec2-deck --test fixtures survey -- --ignored --nocapture`.
#[test]
#[ignore = "development aid"]
fn survey() {
    let mut names: Vec<String> = std::fs::read_dir(fixtures().join("decks"))
        .expect("deck directory")
        .filter_map(|entry| entry.ok()?.file_name().to_str()?.strip_suffix(".nec").map(str::to_string))
        .collect();
    names.sort();
    let filter = std::env::var("NEC2_SURVEY").unwrap_or_default();
    for name in names.iter().filter(|name| name.contains(&filter)) {
        let (deck, reference) = load(name);
        let started = std::time::Instant::now();
        let output = nec2_deck::run(&deck);
        let elapsed = started.elapsed().as_millis();
        let strict = support::compare(&reference, &output.report, None, true);
        let relaxed = support::compare(&reference, &output.report, None, false);
        let status = match &output.error {
            Some(error) => format!("ERROR {error}"),
            None => format!("strict {:4}  relaxed {:4}", strict.len(), relaxed.len()),
        };
        println!("{name:28} {elapsed:5} ms  {status}");
        if !filter.is_empty() {
            let shown = if std::env::var("NEC2_STRICT").is_ok() { &strict } else { &relaxed };
            println!("{}", support::describe(name, shown));
        }
    }
}
