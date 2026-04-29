# Agent Instructions

htdp.mode: autonomous
htdp.transparent: false

## Project conventions

- Use the local `.htdp/prds/` and `.htdp/issues/` artifacts as planning context.
- For the Go v1 triage rewrite, keep the implementation read-only with respect to Nextflow run artifacts.
- Do not invoke `nextflow` from the Go triage core or tests.
