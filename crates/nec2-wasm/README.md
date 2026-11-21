# nec2-wasm

NEC-2, the Numerical Electromagnetics Code for wire antennas (method of moments), in WebAssembly. It is the engine of [nec2-rs](https://github.com/CyberDev-space/nec2-rs), a Rust port of the public-domain nec2c that is checked report by report against it.

You give it a deck of NEC cards. It returns the report in the layout of NEC-2 and the main results as numbers: impedances, currents, the power budget and patterns.

## In the browser

```js
import init, { run } from "nec2-wasm";

await init(); // loads nec2_bg.wasm from the place of the module

const deck = `CM Half-wave dipole in free space
CE
GW 1 21 0 0 -0.25 0 0 0.25 0.001
GE 0
EX 0 1 11 0 1 0
FR 0 1 0 0 299.8 0
RP 0 19 1 1000 0 0 10 0
EN
`;

const output = run(deck);
if (output.error) throw new Error(output.error);

const solution = output.jobs[0].solutions[0];
const [r, x] = solution.sources[0].impedance; // 84.8, 48.0 ohms
const pattern = solution.patterns[0];         // pattern.theta[i], pattern.gainTotal[i]
```

Bundlers that understand `new URL("nec2_bg.wasm", import.meta.url)`, such as Vite, need nothing else.

## In Node

Node cannot fetch a file, so give `init` the bytes:

```js
import { readFile } from "node:fs/promises";
import init, { run } from "nec2-wasm";

const wasm = await readFile(new URL(import.meta.resolve("nec2-wasm/nec2_bg.wasm")));
await init({ module_or_path: wasm });
```

## What it returns

`run(deck)` and `runWith(deck, options)` return an object; the TypeScript types of the package describe every field.

- `report`: the text report.
- `error`: the message of the error that stopped the run, or `null`. An error of the deck does not throw.
- `warnings`: the warnings of the report.
- `jobs`: one entry for every job of the deck (one, unless the deck has NX cards), each with
  - `segments`: tag, ends and radius of every segment, in metres;
  - `solutions`: one for every frequency of a sweep and every direction of an incident wave, with `frequency` (MHz), `sources` and `networkPorts` (voltage, current, impedance, admittance as `[re, im]`, power), `currents` (`[re, im]` in amperes at the centre of every segment), `power` (input, radiated, losses) and `patterns`.

A pattern is stored by columns: `theta[i]`, `phi[i]`, `gainTotal[i]` and so on belong to direction `i`.

The options of `runWith` are `maxEquations`, which refuses a structure that is larger (set it when the deck comes from a user), and `overlap` (`"warn"`, `"error"` or `"off"`), which says what to do with wires that occupy the same space.

The tables of the report for patches, near fields, charges and coupling are in `report` only.

The computation runs in the thread that calls it. In a browser, call it from a Web Worker to keep the page responsive.

## License

MIT. Where the code comes from is described in [PROVENANCE.md](https://github.com/CyberDev-space/nec2-rs/blob/main/PROVENANCE.md).
