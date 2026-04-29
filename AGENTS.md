# Agent Instructions

htdp.mode: autonomous
htdp.transparent: false

## Project conventions

- Use the local `.htdp/prds/` and `.htdp/issues/` artifacts as planning context.
- For the Go v1 triage rewrite, keep the implementation read-only with respect to Nextflow run artifacts.
- Do not invoke `nextflow` from the Go triage core or tests.

## Build and verification

- Build command: `go build ./...`
- Test command: `go test ./...`
- After changing Go source, tests, docs with Go-backed tests, or module files, rebuild and test with `go build ./... && go test ./...` before reporting completion.
