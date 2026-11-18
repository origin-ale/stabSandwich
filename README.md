# Joining CAMPS and Pauli propagation

Simulations of Clifford-augmented MPS (CAMPS) joined with Pauli propagation (PP), the
"sandwich" CAMPS-PP method: a circuit is evolved with CAMPS until the bond dimension reaches
a threshold, and the rest of the circuit is handled by propagating the observable with PP.
Built on [MagicTensors.jl](https://github.com/MagicTensors/MagicTensors.jl),
[PauliPropagation.jl](https://github.com/SparqleSim/PauliPropagation.jl) and
[QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl).

## Layout

| Path | What it is |
|---|---|
| `CampsPP/` | The package: CAMPS evolution and disentangling (submodule `DisentangleCAMPS`), the PauliPropagation ↔ MagicTensors bridge, circuit builders, the `*_circuit_dynamics` drivers, output writers and plain-text progress logging. Self-contained, with its own `Project.toml`/`Manifest.toml`, so it can be shared on its own. |
| `*.jl` (root) | Simulation scripts, see below. |
| `scripts/ScriptSetup.jl` | Shared script setup: `@parameters` (command-line overrides), `param_header`, `setup_threads`. |
| `plotting/` | gnuplot (and one Julia) plot scripts. Most read fetched run folders under `output/`, e.g. `output/carlos_varseed64/`. |
| `scripts/remote/` | Running scripts on the remote machine and fetching the results (`run.sh`, `fetch.sh`; see its README). |
| `test/` | `runtests.jl` (runs the CampsPP suite) and `CONVENTIONS.md` (the conventions the tests pin). |
| `output/` | All results (not tracked by git). |

## Setup

Developed with Julia 1.13.3. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The root environment loads `CampsPP` from the local folder and pins MagicTensors to a
commit (see `[sources]` in `Project.toml`).

## Running a script

```bash
julia --project=. -t 4 Comptime_prob_tm.jl methods=campspp,pp magic_probs=0,0.02 samples=10
```

- **Parameters.** Every value in a script's `@parameters` block can be overridden as
  `key=value`. Values are parsed to the type of the default: numbers (also `pi/4`-style
  expressions; integers also as `1e7`), `true`/`false`, symbols without the colon
  (`camps_crit=entangle`), strings, comma-separated
  lists without spaces (`Ns=12,24`) and ranges (`fractions=0.1:0.1:1.1`). An unknown key is
  an error. `julia --project=. <script>.jl help` lists the parameters and their values.
- **Output.** Each script writes to `output/<script name>/` (override with `outdir=…`). The
  data files written with `initialize_output` start with a `# key=value …` header that
  records all the parameters. Exceptions: `Randrot_complexity.jl`'s `gen_results.dat` has its
  own header (`rng_seed`, `n_samples`, then one block per N), and the per-method
  `comptime_N<N>_<method>_log.txt` trajectories of `Comptime_prob_tm.jl` have none.
- **Progress.** Plain timestamped lines on stdout (no progress bars), so a log redirected
  to a file can be followed with `tail -f`.
- **Threads.** `Comptime_prob_tm.jl` and `TransferredMag_varparams.jl` give BLAS and Strided
  one thread per Julia thread (`-t`). The timing scripts (`Comptime_ratios.jl`,
  `Comptime_3_squarecirc.jl`, `switch_optimization.jl`) fix them at 1.

| Script | Computes |
|---|---|
| `Comptime_prob_tm.jl` | Transferred magnetization in doped XXZ Floquet circuits (N = 46): computation time, resources (bond dimension, number and weight of Paulis), expectation values and switch points vs doping probability, for CAMPS-PP, CAMPS, MPS and PP. |
| `TransferredMag_varparams.jl` | CAMPS-PP transferred-magnetization dynamics over a grid of doping probabilities and domain-wall parameters μ. Written to run samples in parallel threads, but **run it with `-t 1` for now**: with more threads it currently fails ("rotation needs a Hermitian Pauli string"), a thread-safety problem in CampsPP or its dependencies, to be fixed. |
| `switch_optimization.jl` | CAMPS-PP computation time on random rotation circuits as a function of the switch point (how many gates are left to PP). |
| `Comptime_ratios.jl` | Exact CAMPS vs PP computation time on random rotation circuits, vs the aspect ratio t/N at fixed N. |
| `Comptime_3_squarecirc.jl` | Exact CAMPS, PP and MPS computation time on square (t = N) random rotation circuits, vs N. |
| `Randrot_complexity.jl` | Entanglement and stabilizer Rényi entropy after t layers of (random Clifford, random Pauli rotation), vs t/N; also saves two plots. |

Remote runs: `scripts/remote/run.sh -t 10 Comptime_prob_tm.jl samples=25` on the target,
then `scripts/remote/fetch.sh <run id>` locally (see `scripts/remote/README.md`).

## Tests

```bash
julia test/runtests.jl                         # fast tier, ~3 min
STABSANDWICH_SLOW=1 julia test/runtests.jl     # adds the slow checks
```

## References
G. E. Fux, B. Béri, R. Fazio, and E. Tirrito, “Disentangling Magic States with Classically Simulable Quantum Circuits,” Phys. Rev. Lett., vol. 135, no. 26, Dec. 2025, doi: [10.1103/ggp1-byj1](https://doi.org/10.1103/ggp1-byj1).

Z. Liu and B. K. Clark, “Classical simulability of Clifford+T circuits with Clifford-augmented matrix product states,” Aug. 26, 2025, arXiv:2412.17209. doi: [10.48550/arXiv.2412.17209](https://doi.org/10.48550/arXiv.2412.17209).

M. S. Rudolph, T. Jones, Y. Teng, A. Angrisani, and Z. Holmes, "Pauli propagation: A computational framework for simulating quantum systems." Jun. 11, 2025, arXiv:2505.21606. doi: [10.48550/arXiv.2505.21606](https://doi.org/10.48550/arXiv.2505.21606)

E. Rosenberg et al., “Dynamics of magnetization at infinite temperature in a Heisenberg spin chain,” Science, vol. 384, no. 6691, pp. 48–53, Apr. 2024, doi: [10.1126/science.adi7877](https://doi.org/10.1126/science.adi7877).
