# Remote runs

Simulations run on the remote machine, one at a time, and come back through two
SSH jumps. You run and fetch them yourself; Claude sessions can't reach the
remote machine and only read what you fetch into `output/`.

| File | Where it runs | What it does |
|---|---|---|
| `ssh_config.example` | your `~/.ssh/config` | host `stab` through both jumps, with one shared connection |
| `run.sh` | on the target | detached `julia` run with a lock, a run record, and checksums |
| `fetch.sh` | locally | copies a finished run back and verifies it before installing |

## One-time setup

1. Copy the entry from `ssh_config.example` into `~/.ssh/config` and fill in
   the hosts and users.
2. On the target, clone the repo and instantiate it:
   `git clone https://github.com/origin-ale/stabSandwich.git && cd stabSandwich && julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
   If the clone path isn't `~/stabSandwich`, set `STAB_DIR` locally, for example
   `export STAB_DIR=work/stabSandwich`.

## Each working session

```bash
ssh -fN stab                                  # three passwords, once
ssh stab                                      # then, on the target:
  cd stabSandwich
  scripts/remote/run.sh -t 8 Comptime_ratios.jl  # commit and push first; run.sh pulls
  scripts/remote/run.sh --status
  exit
scripts/remote/fetch.sh --list                # back on your machine
scripts/remote/fetch.sh latest                # or a run id from --list
ssh -O exit stab
```

## Run directories

Each run gets `output/runs/<YYYYMMDD-HHMMSS>-<script>/` on the target, and
the same directory locally after a fetch:

- `record.txt`: commit, branch, host, cores, threads, Julia version, command,
  start and end time, exit code, and status (`done`, `failed` or `killed`).
- `log.txt`: Julia's stdout and stderr.
- `files.txt`: the files under `output/` that the run created or modified.
- `files/`: the run's own copy of those files. Scripts write to fixed names, so
  the next run overwrites the copies in `output/`, but not this one.
- `SHA256SUMS`: written last, so a run without it hasn't finished.

`fetch.sh` also copies each output file to its usual path under `output/`, where
the plotting scripts look. If a local file there differs, it keeps the local
file and exits non-zero. Use `-f` to replace it.

## Rules the scripts enforce

- **One run at a time.** A second `run.sh` refuses while a run is alive, which
  keeps timing benchmarks honest and makes it unambiguous which run wrote which
  file. A lock left behind by a dead run is removed automatically.
- **The record describes the code.** `run.sh` refuses to start if tracked files
  have uncommitted changes. It runs `git pull --ff-only` first unless you pass
  `--no-pull`.
- **Nothing unverified enters `output/`.** `fetch.sh` refuses unfinished runs,
  and installs nothing if the transfer is truncated or any checksum fails.
  Re-running it is safe.

To stop a run, `kill` the worker pid that `run.sh --status` shows. Julia is
stopped and the run is recorded as `killed`; its partial outputs are still
listed and fetchable.

Don't write into `output/` on the target by hand during a run: anything
modified in that window is attributed to the run.
