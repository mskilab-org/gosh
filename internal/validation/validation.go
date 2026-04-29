package validation

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/mskilab-org/gosh/internal/cli"
	"github.com/mskilab-org/gosh/internal/domain"
)

type SmokeCommandTiming struct {
	Command string
	Cold    bool
	Elapsed time.Duration
}

type SmokeValidationReport struct {
	RunDir             domain.RunDir
	SyntheticOnly      bool
	FixtureDescription string
	Timings            []SmokeCommandTiming
	Limitations        []string
}

func ValidateSyntheticSmokePerformance(ctx context.Context, runDir domain.RunDir, binaryPath string) (SmokeValidationReport, error) {
	report := SmokeValidationReport{
		SyntheticOnly:      true,
		FixtureDescription: "synthetic trace-backed Nextflow-like fixture; synthetic-only smoke/performance validation, not real nf-gOS validation",
		Limitations: []string{
			"synthetic-only: uses generated trace.tsv and .command.* files instead of a copied/sanitized real nf-gOS run directory",
			"not real nf-gOS validation: trace columns, task volume, filesystem layout, and cache behavior may differ from production runs",
			"uses in-process Execute; does not execute Nextflow, Python, AI, network calls, or an external gosh binary",
			"timings are local smoke measurements, not stable benchmark claims",
		},
	}
	if strings.TrimSpace(binaryPath) != "" {
		report.Limitations = append(report.Limitations, fmt.Sprintf("binaryPath %q was recorded as optional metadata and not executed; in-process Execute was used instead", binaryPath))
	}

	if ctx == nil {
		return report, fmt.Errorf("synthetic smoke validation: nil context")
	}
	if err := ctx.Err(); err != nil {
		return report, fmt.Errorf("synthetic smoke validation: %w", err)
	}

	syntheticPath := strings.TrimSpace(runDir.Path)
	if syntheticPath == "" {
		createdPath, err := os.MkdirTemp("", "gosh-synthetic-smoke-*")
		if err != nil {
			return report, fmt.Errorf("create synthetic smoke run dir: %w", err)
		}
		syntheticPath = createdPath
	} else {
		absolutePath, err := filepath.Abs(syntheticPath)
		if err != nil {
			return report, fmt.Errorf("resolve synthetic smoke run dir %q: %w", runDir.Path, err)
		}
		syntheticPath = absolutePath

		info, err := os.Stat(syntheticPath)
		switch {
		case os.IsNotExist(err):
			if err := os.MkdirAll(syntheticPath, 0o755); err != nil {
				return report, fmt.Errorf("create synthetic smoke run dir %q: %w", syntheticPath, err)
			}
		case err != nil:
			return report, fmt.Errorf("stat synthetic smoke run dir %q: %w", syntheticPath, err)
		case !info.IsDir():
			return report, fmt.Errorf("synthetic smoke run dir %q is not a directory", syntheticPath)
		default:
			entries, err := os.ReadDir(syntheticPath)
			if err != nil {
				return report, fmt.Errorf("read synthetic smoke run dir %q: %w", syntheticPath, err)
			}
			if len(entries) > 0 {
				createdPath, err := os.MkdirTemp(syntheticPath, "gosh-synthetic-smoke-*")
				if err != nil {
					return report, fmt.Errorf("create nested synthetic smoke run dir under %q: %w", syntheticPath, err)
				}
				syntheticPath = createdPath
			}
		}
	}
	syntheticPath = filepath.Clean(syntheticPath)
	report.RunDir = domain.RunDir{Path: syntheticPath}

	if err := ctx.Err(); err != nil {
		return report, fmt.Errorf("synthetic smoke validation: %w", err)
	}

	completedWorkdir := filepath.Join(syntheticPath, "work", "aa", "111111")
	failedWorkdir := filepath.Join(syntheticPath, "work", "bb", "222222")
	cachedWorkdir := filepath.Join(syntheticPath, "work", "cc", "333333")
	for _, workdir := range []string{completedWorkdir, failedWorkdir, cachedWorkdir} {
		if err := os.MkdirAll(workdir, 0o755); err != nil {
			return report, fmt.Errorf("create synthetic workdir %q: %w", workdir, err)
		}
	}

	commandFiles := map[string]string{
		filepath.Join(failedWorkdir, ".command.sh"):  "#!/usr/bin/env bash\necho synthetic failure\nexit 1\n",
		filepath.Join(failedWorkdir, ".command.err"): "synthetic error: not real nf-gOS validation\n",
		filepath.Join(failedWorkdir, ".command.log"): "synthetic log only; no Nextflow execution occurred\nERROR synthetic failure\n",
		filepath.Join(failedWorkdir, ".command.out"): "synthetic stdout\n",
	}
	for path, content := range commandFiles {
		if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
			return report, fmt.Errorf("write synthetic command file %q: %w", path, err)
		}
	}

	traceContent := strings.Join([]string{
		"hash\tstatus\tprocess\tname\ttag\tworkdir\texit\tduration\trealtime\tcpus\tmemory",
		"aa/111111\tCOMPLETED\tSYNTH_ALIGN\tSYNTH_ALIGN (sample-ok)\tsample-ok\t" + completedWorkdir + "\t0\t1m\t60s\t2\t4 GB",
		"bb/222222\tFAILED\tSYNTH_FAIL\tSYNTH_FAIL (sample-failed)\tsample-failed\t" + failedWorkdir + "\t1\t2m\t120s\t4\t8 GB",
		"cc/333333\tCACHED\tSYNTH_CACHE\tSYNTH_CACHE (sample-cached)\tsample-cached\t" + cachedWorkdir + "\t0\t30s\t30s\t1\t1 GB",
		"",
	}, "\n")
	tracePath := filepath.Join(syntheticPath, "trace.tsv")
	if err := os.WriteFile(tracePath, []byte(traceContent), 0o644); err != nil {
		return report, fmt.Errorf("write synthetic trace %q: %w", tracePath, err)
	}

	commands := []struct {
		name string
		cold bool
		args []string
	}{
		{name: "status --json", cold: true, args: []string{"status", "--run-dir", syntheticPath, "--json"}},
		{name: "status --json", cold: false, args: []string{"status", "--run-dir", syntheticPath, "--json"}},
		{name: "tasks --json --status FAILED", cold: false, args: []string{"tasks", "--run-dir", syntheticPath, "--json", "--status", "FAILED"}},
		{name: "inspect bb/222222 --json", cold: false, args: []string{"inspect", "bb/222222", "--run-dir", syntheticPath, "--json"}},
	}

	for _, command := range commands {
		if err := ctx.Err(); err != nil {
			return report, fmt.Errorf("synthetic smoke validation before %q: %w", command.name, err)
		}

		var stdout bytes.Buffer
		var stderr bytes.Buffer
		startedAt := time.Now()
		err := cli.Execute(ctx, command.args, &stdout, &stderr, "synthetic-smoke")
		elapsed := time.Since(startedAt)
		report.Timings = append(report.Timings, SmokeCommandTiming{
			Command: command.name,
			Cold:    command.cold,
			Elapsed: elapsed,
		})
		if err != nil {
			if stderr.Len() > 0 {
				return report, fmt.Errorf("synthetic smoke command %q failed: %w (stderr: %s)", command.name, err, strings.TrimSpace(stderr.String()))
			}
			return report, fmt.Errorf("synthetic smoke command %q failed: %w", command.name, err)
		}
		if strings.TrimSpace(stdout.String()) == "" {
			return report, fmt.Errorf("synthetic smoke command %q produced no output", command.name)
		}
	}

	return report, nil
}
