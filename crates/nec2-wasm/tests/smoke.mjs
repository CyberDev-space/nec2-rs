// Checks the built package in Node: scripts/build-wasm.sh, then
//   node crates/nec2-wasm/tests/smoke.mjs
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const pkg = new URL("../../../target/npm/nec2-wasm/", import.meta.url);
const { default: init, run, runWith, version } = await import(new URL("nec2.js", pkg));
await init({ module_or_path: await readFile(new URL("nec2_bg.wasm", pkg)) });

const manifest = JSON.parse(await readFile(new URL("package.json", pkg), "utf8"));
assert.equal(version(), manifest.version);

const deck = `CM Half-wave dipole in free space
CE
GW 1 21 0 0 -0.25 0 0 0.25 0.001
GE 0
EX 0 1 11 0 1 0
FR 0 3 0 0 290 10
RP 0 19 1 1000 0 0 10 0
EN
`;
const started = performance.now();
const output = run(deck);
const elapsed = performance.now() - started;

assert.equal(output.error, null);
assert.deepEqual(output.warnings, []);
assert.ok(output.report.includes("ANTENNA INPUT PARAMETERS"));
const [job] = output.jobs;
assert.equal(job.segments.length, 21);
assert.deepEqual(job.segments[0].end1, [0, 0, -0.25]);
assert.deepEqual(job.solutions.map((solution) => solution.frequency), [290, 300, 310]);

const solution = job.solutions[1];
const [r, x] = solution.sources[0].impedance;
assert.ok(Math.abs(r - 85.0) < 0.5 && Math.abs(x - 48.5) < 1.0, `impedance ${r} ${x}`);
assert.equal(solution.currents.length, 21);
assert.equal(solution.incidentWave, null);
assert.ok(solution.power.input > 0 && Math.abs(solution.power.radiated - solution.power.input) < 1e-12);
const [pattern] = solution.patterns;
assert.equal(pattern.theta.length, 19);
assert.equal(pattern.theta[9], 90);
assert.ok(Math.abs(pattern.gainTotal[9] - 2.2) < 0.2, `gain ${pattern.gainTotal[9]}`);
assert.equal(pattern.sense[9], "LINEAR");
assert.equal(pattern.averageGain, null);

// An error of the deck is a value, not an exception, and keeps what was computed.
const failed = run(deck.replace("RP 0 19", "RP 7 19"));
assert.match(failed.error, /THE PATTERN MODE MUST BE 0 TO 6/);
assert.equal(failed.jobs[0].solutions.length, 0);

// Options.
const refused = runWith(deck, { maxEquations: 20 });
assert.match(refused.error, /NEEDS 21 EQUATIONS, MORE THAN THE LIMIT OF 20/);
assert.equal(runWith(deck, { overlap: "off" }).error, null);
assert.throws(() => runWith(deck, { overlap: "maybe" }), /overlap must be/);

console.log(`nec2-wasm ${version()}: ok, three frequencies of a 21-segment dipole in ${elapsed.toFixed(1)} ms`);
