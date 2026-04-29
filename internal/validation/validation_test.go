package validation

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
	"github.com/mskilab-org/gosh/internal/run"
)

func TestValidateSyntheticSmokePerformanceRecordsColdAndWarmCLIFlows(t *testing.T) {
	workspace := t.TempDir()

	report, err := ValidateSyntheticSmokePerformance(context.Background(), domain.RunDir{Path: workspace}, "")
	if err != nil {
		t.Fatalf("ValidateSyntheticSmokePerformance returned error: %v", err)
	}

	if !report.SyntheticOnly {
		t.Fatalf("SyntheticOnly = false, want true")
	}
	if report.RunDir.Path != filepath.Clean(workspace) {
		t.Fatalf("RunDir.Path = %q, want %q", report.RunDir.Path, filepath.Clean(workspace))
	}
	for _, want := range []string{"synthetic", "not real nf-gOS"} {
		if !strings.Contains(report.FixtureDescription, want) {
			t.Fatalf("FixtureDescription = %q, want it to contain %q", report.FixtureDescription, want)
		}
	}
	if _, err := os.Stat(run.IndexPath(report.RunDir)); err != nil {
		t.Fatalf("synthetic smoke should build a cold index at %q: %v", run.IndexPath(report.RunDir), err)
	}

	wantCommands := []struct {
		command string
		cold    bool
	}{
		{command: "status --json", cold: true},
		{command: "status --json", cold: false},
		{command: "tasks --json --status FAILED", cold: false},
		{command: "inspect bb/222222 --json", cold: false},
	}
	if len(report.Timings) != len(wantCommands) {
		t.Fatalf("timings length = %d, want %d (%#v)", len(report.Timings), len(wantCommands), report.Timings)
	}
	for index, want := range wantCommands {
		got := report.Timings[index]
		if got.Command != want.command || got.Cold != want.cold {
			t.Fatalf("timing[%d] = %#v, want command %q cold %v", index, got, want.command, want.cold)
		}
		if got.Elapsed < 0 {
			t.Fatalf("timing[%d].Elapsed = %v, want non-negative duration", index, got.Elapsed)
		}
	}

	for _, want := range []string{"synthetic-only", "not real nf-gOS", "Nextflow", "in-process Execute"} {
		if !containsSubstring(report.Limitations, want) {
			t.Fatalf("limitations = %#v, want an entry containing %q", report.Limitations, want)
		}
	}
}

func TestValidateSyntheticSmokePerformanceTreatsBinaryPathAsOptionalMetadata(t *testing.T) {
	workspace := t.TempDir()
	missingBinary := filepath.Join(workspace, "missing-gosh-binary")

	report, err := ValidateSyntheticSmokePerformance(context.Background(), domain.RunDir{Path: workspace}, missingBinary)
	if err != nil {
		t.Fatalf("ValidateSyntheticSmokePerformance with missing binaryPath returned error: %v", err)
	}
	if len(report.Timings) == 0 {
		t.Fatalf("timings length = 0, want synthetic CLI commands to run in-process")
	}
	for _, want := range []string{missingBinary, "not executed"} {
		if !containsSubstring(report.Limitations, want) {
			t.Fatalf("limitations = %#v, want an entry containing %q", report.Limitations, want)
		}
	}
}

func TestValidateSyntheticSmokePerformanceHonorsCanceledContext(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	report, err := ValidateSyntheticSmokePerformance(ctx, domain.RunDir{Path: t.TempDir()}, "")
	if err == nil {
		t.Fatalf("ValidateSyntheticSmokePerformance with canceled context returned nil error")
	}
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("error = %v, want context.Canceled", err)
	}
	if !report.SyntheticOnly {
		t.Fatalf("SyntheticOnly = false, want partial report to still identify synthetic-only validation")
	}
	if len(report.Timings) != 0 {
		t.Fatalf("timings = %#v, want no commands after pre-canceled context", report.Timings)
	}
}

func containsSubstring(values []string, want string) bool {
	for _, value := range values {
		if strings.Contains(value, want) {
			return true
		}
	}
	return false
}
