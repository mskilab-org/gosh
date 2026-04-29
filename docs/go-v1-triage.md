# Go v1 fast triage CLI

The Go v1 CLI is a focused, read-only triage surface for existing Nextflow run artifacts. It is not a feature-preserving port of the old Python CLI. The authoritative v1 commands are `status`, `tasks`, `inspect`, and `index`.

Use v1 when you have an already-created run directory and want to answer:

- what failed,
- which canonical task/workdir ID identifies it,
- which trace-backed task rows match a filter, and
- which `.command.*` files are worth reading next.

## Scope and safety boundary

V1 is deliberately narrow:

- It reads selected trace files, selected Nextflow logs, and bounded `.command.*` snippets from task workdirs.
- It may create or update the run-local cache at `.gosh/index.sqlite`.
- It does not mutate trace files, `.nextflow.log` files, task workdirs, `.command.*` files, pipeline outputs, samplesheets, params files, or Nextflow run metadata.
- It does not execute Nextflow, does not invoke Python, does not call AI services, and does not use the network.
- `inspect` inventories command files and reads snippets; it never executes `.command.sh` or `.command.run`.

This means the CLI is read-only with respect to Nextflow run artifacts, with one explicit local write boundary: the `.gosh/index.sqlite` cache used for fast repeated triage.

## Command surface

```text
gosh status  [--run-dir DIR] [--json]
gosh tasks   [--run-dir DIR] [filters] [--json]
gosh inspect <task> [--run-dir DIR] [--json]
gosh index   [--run-dir DIR] [--refresh] [--json]
```

Global options can appear before the command or as command-local options:

- `--run-dir` / `-d`: select the Nextflow run directory; default is `.`.
- `--json`: emit JSON instead of human text.
- `--format human|json`: explicit output format.

Examples using the long run-dir flag:

```sh
gosh status --run-dir /path/to/run
gosh tasks --run-dir /path/to/run --status FAILED
gosh inspect ab/c123def --run-dir /path/to/run
```

Examples using `-d` and JSON output:

```sh
gosh status -d /path/to/run --json
gosh tasks -d /path/to/run --json --status FAILED
gosh inspect ab/c123def -d /path/to/run --json
```

## Trace-first artifact discovery

For each selected run directory, v1 discovers supported artifacts in this order:

1. Trace candidates: `trace*.txt`, `trace*.csv`, `trace*.tsv`.
2. Log candidates: `.nextflow.log`, `.nextflow_*.log`.

If any trace candidate exists, v1 selects the newest trace file as the primary source and enters `trace-backed` mode. If a log is also present, the newest log is paired as a secondary source for metadata/diagnostics.

If no trace exists but a supported log exists, v1 enters `log-only` mode. If neither a trace nor a supported log exists, commands report `unsupported_artifacts` diagnostics and recommend future runs use `-with-trace`.

## Cache and `gosh index`

The run-local cache path is `.gosh/index.sqlite` under the selected run directory. Cache metadata records the selected source paths, mtimes, sizes, mode, schema version, build time, freshness, and task count.

Trace-backed `status`, `tasks`, and `inspect` call the fresh-index path automatically. When the selected trace/log fingerprints are missing or stale, those commands rebuild the trace-backed task rows before rendering output.

`gosh index` is a diagnostics command:

```sh
gosh index --run-dir /path/to/run
gosh index --run-dir /path/to/run --refresh
gosh index --refresh
```

Current behavior details:

- `gosh index` without `--refresh` reports diagnostics and does not create the cache just to inspect it.
- `gosh index --run-dir /path/to/run --refresh` creates or updates `.gosh/index.sqlite` metadata for the selected artifact set.
- `gosh index --refresh` is metadata-oriented in the current implementation; it does not parse task rows the way `status`, `tasks`, and `inspect` do when they rebuild a trace-backed task index.
- In unsupported mode, `gosh index --refresh` can still write unsupported metadata so diagnostics are explicit and repeatable.

## Canonical IDs and selectors

The canonical task ID is the Nextflow workdir hash path, for example `ab/c123def`. V1 derives this ID from a trace `hash` column, a trace `workdir` column, or a full workdir path when possible. IDs are normalized to lower-case `xx/rest` form.

`gosh tasks` always prints the canonical ID as the first column. `gosh inspect` accepts:

- a canonical ID such as `ab/c123def`,
- a full workdir path, or
- a human selector that matches exactly one task by process/name/tag.

Ambiguous human selectors render a disambiguation table instead of guessing. Missing selectors render a not-found diagnostic. V1 does not fabricate task rows or command-file paths when the trace does not provide enough information.

## Trace-backed workflow: `status` -> `tasks` -> `inspect`

Start with `status` to discover the mode, selected sources, freshness, counts, and failed preview:

```sh
gosh status --run-dir /path/to/run
```

In trace-backed mode, human output includes fields such as:

```text
run_dir: /path/to/run
mode: trace-backed
index: /path/to/run/.gosh/index.sqlite
freshness: fresh
counts:
  COMPLETED: 10
  FAILED: 1
failed_count: 1
failed_preview:
  - id=ab/c123def status=FAILED process=ALIGN name=ALIGN (sample-01) tag=sample-01 workdir=/path/to/run/work/ab/c123def exit=137 error=-
```

Use `tasks` to list rows or narrow to the failed row(s):

```sh
gosh tasks --run-dir /path/to/run --status FAILED
gosh tasks --run-dir /path/to/run --process ALIGN
gosh tasks --run-dir /path/to/run --sample sample-01
```

Task filters are case-insensitive where applicable. `--process` and `--name` are substring filters. `--status` uses normalized exact status matching when possible. `--sample` is a name/tag substring filter, not a parsed biological sample field.

The default task table is source-order and starts with:

```text
id	status	process	name/tag	workdir	exit	duration	realtime	cpus	memory
```

Use the canonical ID from `status` or `tasks` with `inspect`:

```sh
gosh inspect ab/c123def --run-dir /path/to/run
```

For an exact trace-backed match, `inspect` renders task metadata plus command-file inventory and bounded snippets for expected files such as `.command.sh`, `.command.log`, `.command.err`, `.command.out`, and `.command.run` when they exist. Missing command files are reported as missing; they are not created.

## JSON examples

All v1 commands support JSON with `--json` or `--format json`:

```sh
gosh status -d /path/to/run --json
gosh tasks -d /path/to/run --json --status FAILED
gosh inspect ab/c123def -d /path/to/run --json
gosh index -d /path/to/run --json
```

JSON output is intended for automation and golden tests. It includes stable top-level markers such as `"format": "json"`, selected mode (`"trace-backed"`, `"log-only"`, or `"unsupported"`), query metadata, diagnostics, and command-specific payloads.

## Log-only degraded behavior

A log-only run is a run where no supported trace file exists but a supported Nextflow log exists. In that mode, v1 is intentionally conservative.

`gosh status` parses only deterministic failure evidence from the selected log. It can show parseable failures, workdirs, process/name text, exit code, and a short error summary when those appear in recognizable Nextflow error blocks. It also states that complete counts are unavailable:

```text
mode: log-only
counts: unavailable (log-only mode; trace file required)
failed_count: 1
diagnostics:
  - warning log_only_degraded: complete task/resource/status data is unavailable without a Nextflow trace file
  - info nextflow_with_trace_recommended: Run future Nextflow workflows with -with-trace
```

If a log exists but no parseable process failure is found, `status` still reports `mode: log-only`, `failed_count: 0`, and a `log_only_no_parseable_failures` diagnostic. This means “no deterministic failure block was parsed,” not “the workflow succeeded.”

`gosh tasks` and `gosh inspect` require a trace-backed task index. In log-only mode they return diagnostics such as `tasks_unavailable_log_only` or `inspect_unavailable_log_only` and do not fabricate task tables, selectors, dossiers, or command-file paths from incomplete log evidence.

For future runs, prefer:

```sh
nextflow run ... -with-trace
```

The resulting trace file lets v1 build a complete task index.

## Unsupported artifact behavior

If the selected run directory has neither a supported trace nor a supported Nextflow log, v1 reports searched patterns and `unsupported_artifacts`. `status`, `tasks`, and `inspect` return command errors because there is no usable triage source. `index` can still report diagnostics, and `index --refresh` can persist unsupported metadata for repeatable diagnostics.

## Non-compatibility with old Python commands

The Go v1 command surface is intentionally not compatible with old Python command groups. `gosh run`, `gosh debug`, and `gosh help` are not v1 commands. Use `gosh --help` and `gosh <command> --help` for v1 help.

Old behaviors such as pipeline launch/orchestration, module loading, parameter mutation, run-output management, `gosh debug log`, AI explanation/help, and Python-backed Nextflow shell-outs are outside v1. During the branch phase, the Python implementation may remain in the repository as a separate/legacy concern, but the Go triage core does not depend on it.

## Synthetic-only validation limitation

`ValidateSyntheticSmokePerformance` validates only a synthetic CLI smoke path: cold `status --json`, warm `status --json`, `tasks --json --status FAILED`, and `inspect bb/222222 --json` against generated trace/workdir fixtures. It uses in-process `Execute` and generated `.command.*` files.

This is synthetic-only, not real nf-gOS validation:

- No sanitized real nf-gOS run fixture was provided.
- The synthetic trace columns, task volume, filesystem layout, and cache behavior may differ from production nf-gOS runs.
- The validation does not execute Nextflow, does not invoke Python, does not call AI services, does not use the network, and does not execute an external `gosh` binary.
- Recorded timings are local smoke measurements, not benchmark claims.

Therefore these docs do not claim real nf-gOS performance numbers. The qualitative design goal is fast repeated trace-backed triage from a local SQLite cache, but real-run performance and parser coverage still require a copied/sanitized nf-gOS run fixture.
