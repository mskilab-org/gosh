package nflog

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestParseLogOnlyFailuresParsesAndNormalizesSelectedLog(t *testing.T) {
	runPath := t.TempDir()
	workdir := filepath.Join(runPath, "work", "ab", "c123def")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", workdir, err)
	}

	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:01:00.000 [Task monitor] ERROR nextflow.processor.TaskProcessor - Error executing process > 'PIPE:ALIGN (sample-01)'",
		"",
		"Command exit status:",
		"  137",
		"",
		"Command error:",
		"  killed by scheduler",
		"",
		"Work dir:",
		"  AB/C123DEF",
		"Apr-28 12:01:01.000 [main] INFO nextflow.Session - cleanup after first failure",
		"ERROR ~ Error executing process > 'PIPE:QC'",
		"",
		"Command error:",
		"  task failed before workdir was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyFailures(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyFailures returned error: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("ParseLogOnlyFailures returned %d failures, want 2: %#v", len(got), got)
	}

	first := got[0]
	if first.ID != "ab/c123def" || first.Workdir != filepath.Clean(workdir) {
		t.Fatalf("first ID/Workdir = %q/%q, want ab/c123def/%q", first.ID, first.Workdir, filepath.Clean(workdir))
	}
	if first.Process != "PIPE:ALIGN" || first.Name != "sample-01" {
		t.Fatalf("first Process/Name = %q/%q, want PIPE:ALIGN/sample-01", first.Process, first.Name)
	}
	if first.Exit == nil || *first.Exit != 137 {
		if first.Exit == nil {
			t.Fatalf("first Exit = nil, want 137")
		}
		t.Fatalf("first Exit = %d, want 137", *first.Exit)
	}
	if first.ErrorSummary != "killed by scheduler" {
		t.Fatalf("first ErrorSummary = %q, want deterministic command error summary", first.ErrorSummary)
	}
	if !strings.Contains(first.ErrorBlock, "Error executing process > 'PIPE:ALIGN (sample-01)'") {
		t.Fatalf("first ErrorBlock = %q, want original failure block", first.ErrorBlock)
	}

	second := got[1]
	if second.ID != "" || second.Workdir != "" || second.Exit != nil {
		t.Fatalf("second ID/Workdir/Exit = %q/%q/%#v, want partial evidence with empty workdir and nil exit", second.ID, second.Workdir, second.Exit)
	}
	if second.Process != "PIPE:QC" || second.Name != "" {
		t.Fatalf("second Process/Name = %q/%q, want PIPE:QC/empty", second.Process, second.Name)
	}
	if second.ErrorSummary != "task failed before workdir was reported" {
		t.Fatalf("second ErrorSummary = %q, want deterministic command error summary", second.ErrorSummary)
	}
}

func TestParseLogOnlyFailuresReturnsEmptySliceWhenNoFailureEvidenceAppears(t *testing.T) {
	runPath := t.TempDir()
	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:02:00.000 [main] ERROR nextflow.Session - Pipeline aborted before any process failure was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyFailures(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyFailures returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ParseLogOnlyFailures returned nil, want empty slice for no log-only evidence")
	}
	if len(got) != 0 {
		t.Fatalf("ParseLogOnlyFailures returned %#v, want no log-only failures", got)
	}
}

func TestParseLogOnlyFailuresRejectsInvalidInputs(t *testing.T) {
	runPath := t.TempDir()
	logPath := filepath.Join(runPath, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("ERROR ~ Error executing process > 'PIPE:QC'\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	t.Run("nil context", func(t *testing.T) {
		got, err := ParseLogOnlyFailures(nil, domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyFailures(nil context) returned nil error and failures %#v", got)
		}
		if got != nil {
			t.Fatalf("failures on error = %#v, want nil", got)
		}
		if !strings.Contains(err.Error(), "nil context") {
			t.Fatalf("error = %q, want it to mention nil context", err.Error())
		}
	})

	t.Run("canceled context", func(t *testing.T) {
		ctx, cancel := context.WithCancel(context.Background())
		cancel()

		got, err := ParseLogOnlyFailures(ctx, domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyFailures(canceled context) returned nil error and failures %#v", got)
		}
		if got != nil {
			t.Fatalf("failures on error = %#v, want nil", got)
		}
		if !errors.Is(err, context.Canceled) {
			t.Fatalf("error = %v, want context.Canceled", err)
		}
	})

	t.Run("wrong source kind", func(t *testing.T) {
		missingPath := filepath.Join(runPath, "missing-trace.txt")
		got, err := ParseLogOnlyFailures(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: missingPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyFailures(wrong source kind) returned nil error and failures %#v", got)
		}
		if got != nil {
			t.Fatalf("failures on error = %#v, want nil", got)
		}
		for _, want := range []string{"parse log-only failures", "log source", string(domain.SourceKindTrace), string(domain.SourceKindLog)} {
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
			}
		}
	})

	t.Run("missing selected log", func(t *testing.T) {
		missingPath := filepath.Join(runPath, "missing.nextflow.log")
		got, err := ParseLogOnlyFailures(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: missingPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyFailures(missing selected log) returned nil error and failures %#v", got)
		}
		if got != nil {
			t.Fatalf("failures on error = %#v, want nil", got)
		}
		for _, want := range []string{"parse log-only failures", "open source", missingPath} {
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
			}
		}
	})
}

func TestParseLogOnlyTaskEvidenceCombinesLifecycleFailureAndWorkdirEvidence(t *testing.T) {
	runPath := t.TempDir()
	failedWorkdir := filepath.Join(runPath, "work", "ab", "c123def")
	if err := os.MkdirAll(failedWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", failedWorkdir, err)
	}
	if err := os.WriteFile(filepath.Join(failedWorkdir, ".exitcode"), []byte("137\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(.exitcode) returned error: %v", err)
	}
	commandErr := "stderr evidence from referenced workdir"
	if err := os.WriteFile(filepath.Join(failedWorkdir, string(domain.CommandFileErr)), []byte(commandErr), 0o644); err != nil {
		t.Fatalf("WriteFile(.command.err) returned error: %v", err)
	}

	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:00.000 [Task submitter] INFO nextflow.Session - [AB/C123DEF] Submitted process > PIPE:ALIGN (sample-01)",
		"Apr-28 12:00:30.000 [Task monitor] DEBUG nextflow.processor.TaskPollingMonitor - [CD/EF456] Completed process > PIPE:QC (sample-02)",
		"Apr-28 12:01:00.000 [Task monitor] ERROR nextflow.processor.TaskProcessor - Error executing process > 'PIPE:ALIGN (sample-01)'",
		"",
		"Command exit status:",
		"  137",
		"",
		"Command error:",
		"  killed by scheduler",
		"",
		"Work dir:",
		"  AB/C123DEF",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyTaskEvidence returned error: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("ParseLogOnlyTaskEvidence returned %d rows, want combined failed row plus completed row: %#v", len(got), got)
	}

	failed := got[0]
	if failed.ID != "ab/c123def" || failed.Workdir != filepath.Clean(failedWorkdir) {
		t.Fatalf("failed ID/Workdir = %q/%q, want ab/c123def/%q", failed.ID, failed.Workdir, filepath.Clean(failedWorkdir))
	}
	if failed.Process != "PIPE:ALIGN" || failed.Name != "sample-01" || failed.ObservedStatus != domain.TaskStatusFailed {
		t.Fatalf("failed Process/Name/Status = %q/%q/%q, want PIPE:ALIGN/sample-01/FAILED", failed.Process, failed.Name, failed.ObservedStatus)
	}
	if failed.Exit == nil || *failed.Exit != 137 {
		if failed.Exit == nil {
			t.Fatalf("failed Exit = nil, want 137")
		}
		t.Fatalf("failed Exit = %d, want 137", *failed.Exit)
	}
	if failed.ErrorSummary != "killed by scheduler" || !strings.Contains(failed.ErrorBlock, "Error executing process > 'PIPE:ALIGN (sample-01)'") {
		t.Fatalf("failed error fields = summary %q block %q, want final failure-block evidence preserved", failed.ErrorSummary, failed.ErrorBlock)
	}
	if !failed.CommandFilesAvailable || failed.Completeness != domain.LogOnlyEvidencePartial {
		t.Fatalf("failed command availability/completeness = %v/%q, want true/partial after workdir enrichment", failed.CommandFilesAvailable, failed.Completeness)
	}
	for _, want := range []domain.LogOnlyEvidenceSourceKind{domain.LogOnlyEvidenceSourceLog, domain.LogOnlyEvidenceSourceWorkdir, domain.LogOnlyEvidenceSourceCommand} {
		found := false
		for _, source := range failed.Sources {
			if source.Kind == want {
				found = true
				break
			}
		}
		if !found {
			t.Fatalf("failed Sources = %#v, want source kind %q", failed.Sources, want)
		}
	}

	completed := got[1]
	if completed.ID != "cd/ef456" || completed.Process != "PIPE:QC" || completed.Name != "sample-02" || completed.ObservedStatus != domain.TaskStatusCompleted {
		t.Fatalf("completed evidence = %#v, want lifecycle-only completed PIPE:QC sample-02", completed)
	}
	if completed.Workdir != "" || completed.Exit != nil || completed.CommandFilesAvailable || completed.ErrorSummary != "" {
		t.Fatalf("completed evidence = %#v, want no fabricated workdir/exit/error details", completed)
	}
}

func TestParseLogOnlyTaskEvidenceEnrichesLifecycleHashWorkdir(t *testing.T) {
	runPath := t.TempDir()
	workdir := filepath.Join(runPath, "work", "aa", "111111")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", workdir, err)
	}
	if err := os.WriteFile(filepath.Join(workdir, ".exitcode"), []byte("2\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(.exitcode) returned error: %v", err)
	}
	commandErr := "fatal error from lifecycle-only workdir"
	if err := os.WriteFile(filepath.Join(workdir, string(domain.CommandFileErr)), []byte(commandErr), 0o644); err != nil {
		t.Fatalf("WriteFile(.command.err) returned error: %v", err)
	}
	unreferencedWorkdir := filepath.Join(runPath, "work", "bb", "222222")
	if err := os.MkdirAll(unreferencedWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", unreferencedWorkdir, err)
	}
	if err := os.WriteFile(filepath.Join(unreferencedWorkdir, ".exitcode"), []byte("99\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(unreferenced .exitcode) returned error: %v", err)
	}

	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:00.000 [Task submitter] INFO nextflow.Session - [AA/111111] Submitted process > PIPE:ONLY (sample-a)",
		"Apr-28 12:00:30.000 [Task monitor] ERROR nextflow.processor.TaskPollingMonitor - [aa/111111] Failed process > PIPE:ONLY (sample-a)",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyTaskEvidence returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("ParseLogOnlyTaskEvidence returned %d rows, want one lifecycle evidence row: %#v", len(got), got)
	}

	evidence := got[0]
	if evidence.ID != "aa/111111" || evidence.Workdir != filepath.Clean(workdir) {
		t.Fatalf("ID/Workdir = %q/%q, want aa/111111/%q", evidence.ID, evidence.Workdir, filepath.Clean(workdir))
	}
	if evidence.ObservedStatus != domain.TaskStatusFailed || evidence.Process != "PIPE:ONLY" || evidence.Name != "sample-a" {
		t.Fatalf("status/process/name = %q/%q/%q, want FAILED/PIPE:ONLY/sample-a", evidence.ObservedStatus, evidence.Process, evidence.Name)
	}
	if evidence.Exit == nil || *evidence.Exit != 2 {
		if evidence.Exit == nil {
			t.Fatalf("Exit = nil, want 2 from referenced .exitcode")
		}
		t.Fatalf("Exit = %d, want 2 from referenced .exitcode", *evidence.Exit)
	}
	if evidence.ErrorBlock != commandErr || evidence.ErrorSummary != commandErr {
		t.Fatalf("ErrorBlock/ErrorSummary = %q/%q, want bounded command error %q", evidence.ErrorBlock, evidence.ErrorSummary, commandErr)
	}
	if !evidence.CommandFilesAvailable {
		t.Fatalf("CommandFilesAvailable = false, want true for referenced workdir command evidence")
	}
	for _, source := range evidence.Sources {
		if source.Path == filepath.Join(unreferencedWorkdir, ".exitcode") {
			t.Fatalf("Sources = %#v, should not include unreferenced workdir evidence", evidence.Sources)
		}
	}
}

func TestParseLogOnlyTaskEvidenceResolvesShortLifecycleHashToLongWorkdir(t *testing.T) {
	runPath := t.TempDir()
	workdir := filepath.Join(runPath, "work", "9e", "c300c50150c213c2d44ca9e4624d8c")
	if err := os.MkdirAll(workdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", workdir, err)
	}
	if err := os.WriteFile(filepath.Join(workdir, ".exitcode"), []byte("137\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(.exitcode) returned error: %v", err)
	}
	commandErr := "fatal error from prefix-resolved lifecycle workdir"
	if err := os.WriteFile(filepath.Join(workdir, string(domain.CommandFileErr)), []byte(commandErr), 0o644); err != nil {
		t.Fatalf("WriteFile(.command.err) returned error: %v", err)
	}

	// A same-suffix directory in a different shard must not be considered; only work/9e is relevant.
	unreferencedWorkdir := filepath.Join(runPath, "work", "00", "c300c50150c213c2d44ca9e4624d8c")
	if err := os.MkdirAll(unreferencedWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", unreferencedWorkdir, err)
	}
	if err := os.WriteFile(filepath.Join(unreferencedWorkdir, ".exitcode"), []byte("99\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(unreferenced .exitcode) returned error: %v", err)
	}

	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:30.000 [Task monitor] ERROR nextflow.processor.TaskPollingMonitor - [9e/c300c5] failed process > PIPE:STEP (sample)",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyTaskEvidence returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("ParseLogOnlyTaskEvidence returned %d rows, want one lifecycle evidence row: %#v", len(got), got)
	}

	evidence := got[0]
	if evidence.ID != "9e/c300c5" || evidence.Workdir != filepath.Clean(workdir) {
		t.Fatalf("ID/Workdir = %q/%q, want 9e/c300c5/%q", evidence.ID, evidence.Workdir, filepath.Clean(workdir))
	}
	if evidence.ObservedStatus != domain.TaskStatusFailed || evidence.Process != "PIPE:STEP" || evidence.Name != "sample" {
		t.Fatalf("status/process/name = %q/%q/%q, want FAILED/PIPE:STEP/sample", evidence.ObservedStatus, evidence.Process, evidence.Name)
	}
	if evidence.Exit == nil || *evidence.Exit != 137 {
		if evidence.Exit == nil {
			t.Fatalf("Exit = nil, want 137 from prefix-resolved .exitcode")
		}
		t.Fatalf("Exit = %d, want 137 from prefix-resolved .exitcode", *evidence.Exit)
	}
	if evidence.ErrorBlock != commandErr || evidence.ErrorSummary != commandErr {
		t.Fatalf("ErrorBlock/ErrorSummary = %q/%q, want bounded command error %q", evidence.ErrorBlock, evidence.ErrorSummary, commandErr)
	}
	if !evidence.CommandFilesAvailable {
		t.Fatalf("CommandFilesAvailable = false, want true for prefix-resolved workdir command evidence")
	}
	for _, wantPath := range []string{filepath.Clean(workdir), filepath.Join(workdir, ".exitcode"), filepath.Join(workdir, string(domain.CommandFileErr))} {
		found := false
		for _, source := range evidence.Sources {
			if source.Path == wantPath {
				found = true
				break
			}
		}
		if !found {
			t.Fatalf("Sources = %#v, want source path %q", evidence.Sources, wantPath)
		}
	}
	for _, source := range evidence.Sources {
		if strings.Contains(source.Path, unreferencedWorkdir) {
			t.Fatalf("Sources = %#v, should not include unreferenced shard workdir %q", evidence.Sources, unreferencedWorkdir)
		}
	}
}

func TestParseLogOnlyTaskEvidenceKeepsFinalFailureOnlyEvidence(t *testing.T) {
	runPath := t.TempDir()
	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"ERROR ~ Error executing process > 'PIPE:QC'",
		"",
		"Command error:",
		"  task failed before workdir was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyTaskEvidence returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("ParseLogOnlyTaskEvidence returned %d rows, want one partial final-failure row: %#v", len(got), got)
	}
	if got[0].ID != "" || got[0].Workdir != "" || got[0].Exit != nil {
		t.Fatalf("ID/Workdir/Exit = %q/%q/%#v, want missing fields preserved for partial final-failure evidence", got[0].ID, got[0].Workdir, got[0].Exit)
	}
	if got[0].Process != "PIPE:QC" || got[0].Name != "" || got[0].ObservedStatus != domain.TaskStatusFailed {
		t.Fatalf("Process/Name/Status = %q/%q/%q, want PIPE:QC/empty/FAILED", got[0].Process, got[0].Name, got[0].ObservedStatus)
	}
	if got[0].ErrorSummary != "task failed before workdir was reported" {
		t.Fatalf("ErrorSummary = %q, want deterministic command error summary", got[0].ErrorSummary)
	}
}

func TestParseLogOnlyTaskEvidenceReturnsEmptySliceWhenNoEvidenceAppears(t *testing.T) {
	runPath := t.TempDir()
	logPath := filepath.Join(runPath, ".nextflow.log")
	content := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:02:00.000 [main] ERROR nextflow.Session - Pipeline aborted before task evidence was reported",
		"",
	}, "\n")
	if err := os.WriteFile(logPath, []byte(content), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
	if err != nil {
		t.Fatalf("ParseLogOnlyTaskEvidence returned error: %v", err)
	}
	if got == nil {
		t.Fatalf("ParseLogOnlyTaskEvidence returned nil, want empty slice for no observed evidence")
	}
	if len(got) != 0 {
		t.Fatalf("ParseLogOnlyTaskEvidence returned %#v, want no evidence", got)
	}
}

func TestParseLogOnlyTaskEvidenceRejectsInvalidInputs(t *testing.T) {
	runPath := t.TempDir()
	logPath := filepath.Join(runPath, ".nextflow.log")
	if err := os.WriteFile(logPath, []byte("[ab/c123def] Failed process > PIPE:QC\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(%q) returned error: %v", logPath, err)
	}

	t.Run("nil context", func(t *testing.T) {
		got, err := ParseLogOnlyTaskEvidence(nil, domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyTaskEvidence(nil context) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on error = %#v, want nil", got)
		}
		if !strings.Contains(err.Error(), "nil context") {
			t.Fatalf("error = %q, want it to mention nil context", err.Error())
		}
	})

	t.Run("canceled context", func(t *testing.T) {
		ctx, cancel := context.WithCancel(context.Background())
		cancel()

		got, err := ParseLogOnlyTaskEvidence(ctx, domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: logPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyTaskEvidence(canceled context) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on error = %#v, want nil", got)
		}
		if !errors.Is(err, context.Canceled) {
			t.Fatalf("error = %v, want context.Canceled", err)
		}
	})

	t.Run("wrong source kind", func(t *testing.T) {
		got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindTrace, Path: logPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyTaskEvidence(wrong source kind) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on error = %#v, want nil", got)
		}
		for _, want := range []string{"parse log-only task evidence", "log source", string(domain.SourceKindTrace), string(domain.SourceKindLog)} {
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
			}
		}
	})

	t.Run("missing selected log", func(t *testing.T) {
		missingPath := filepath.Join(runPath, "missing.nextflow.log")
		got, err := ParseLogOnlyTaskEvidence(context.Background(), domain.RunDir{Path: runPath}, domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: missingPath})
		if err == nil {
			t.Fatalf("ParseLogOnlyTaskEvidence(missing selected log) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on error = %#v, want nil", got)
		}
		for _, want := range []string{"parse log-only task evidence", "open source", missingPath} {
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
			}
		}
	})
}

func TestExtractLifecycleEvidenceParsesHashPrefixedLifecycleLines(t *testing.T) {
	input := strings.Join([]string{
		"Apr-28 12:00:00.000 [Task submitter] INFO nextflow.Session - [AB/C123DEF] Submitted process > PIPE:ALIGN (sample-01)",
		"Apr-28 12:01:00.000 [Task monitor] DEBUG nextflow.processor.TaskPollingMonitor - [ab/c123def] Completed process > PIPE:ALIGN (sample-01)",
		"Apr-28 12:02:00.000 [Task submitter] INFO nextflow.Session - [DE/F456] Cached process > PIPE:CACHE (sample-cached)",
		"Apr-28 12:03:00.000 [Task monitor] ERROR nextflow.processor.TaskPollingMonitor - [12/ABCDEF] Failed process > PIPE:QUANT (tumor-02)",
		"Apr-28 12:04:00.000 [Task submitter] INFO nextflow.Session - [34/BBBBBB] Submitted process > PIPE:WAITING",
		"Apr-28 12:05:00.000 [main] ERROR nextflow.Session - Pipeline aborted after task failure",
		"",
	}, "\n")

	got, err := ExtractLifecycleEvidence(strings.NewReader(input))
	if err != nil {
		t.Fatalf("ExtractLifecycleEvidence returned error: %v", err)
	}
	if len(got) != 4 {
		t.Fatalf("ExtractLifecycleEvidence returned %d evidence rows, want 4: %#v", len(got), got)
	}

	first := got[0]
	if first.ID != "ab/c123def" || first.Process != "PIPE:ALIGN" || first.Name != "sample-01" || first.ObservedStatus != domain.TaskStatusCompleted {
		t.Fatalf("first evidence = %#v, want completed PIPE:ALIGN sample evidence with canonical id ab/c123def", first)
	}
	if first.Workdir != "" || first.Exit != nil || first.ErrorSummary != "" || first.ErrorBlock != "" {
		t.Fatalf("first evidence has fabricated fields: %#v, want no workdir/exit/error details from lifecycle lines", first)
	}
	if first.Completeness != domain.LogOnlyEvidencePartial || first.CommandFilesAvailable {
		t.Fatalf("first completeness/command availability = %q/%v, want partial/false", first.Completeness, first.CommandFilesAvailable)
	}
	if len(first.Sources) != 2 {
		t.Fatalf("first Sources = %#v, want submitted and completed log source markers", first.Sources)
	}
	for _, source := range first.Sources {
		if source.Kind != domain.LogOnlyEvidenceSourceLog {
			t.Fatalf("first source kind = %q, want log source marker", source.Kind)
		}
	}

	second := got[1]
	if second.ID != "de/f456" || second.Process != "PIPE:CACHE" || second.Name != "sample-cached" || second.ObservedStatus != domain.TaskStatusCached {
		t.Fatalf("second evidence = %#v, want cached PIPE:CACHE sample evidence", second)
	}

	third := got[2]
	if third.ID != "12/abcdef" || third.Process != "PIPE:QUANT" || third.Name != "tumor-02" || third.ObservedStatus != domain.TaskStatusFailed {
		t.Fatalf("third evidence = %#v, want failed PIPE:QUANT tumor evidence", third)
	}

	fourth := got[3]
	if fourth.ID != "34/bbbbbb" || fourth.Process != "PIPE:WAITING" || fourth.Name != "" || fourth.ObservedStatus != domain.TaskStatusSubmitted {
		t.Fatalf("fourth evidence = %#v, want submitted untagged PIPE:WAITING evidence", fourth)
	}
}

func TestExtractLifecycleEvidenceReturnsEmptySliceWhenNoLifecycleEvidenceAppears(t *testing.T) {
	inputs := []string{
		"",
		"   \r\n\t\r\n",
		strings.Join([]string{
			"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
			"ERROR ~ Error executing process > 'PIPE:QC'",
			"[not/a-hash] Submitted process > PIPE:BOGUS",
			"[ab/c123def] process > PIPE:PROGRESS [100%] 1 of 1",
		}, "\n"),
	}

	for _, input := range inputs {
		t.Run(input, func(t *testing.T) {
			got, err := ExtractLifecycleEvidence(strings.NewReader(input))
			if err != nil {
				t.Fatalf("ExtractLifecycleEvidence(%q) returned error: %v", input, err)
			}
			if got == nil {
				t.Fatalf("ExtractLifecycleEvidence(%q) returned nil, want empty slice", input)
			}
			if len(got) != 0 {
				t.Fatalf("ExtractLifecycleEvidence(%q) returned %#v, want no lifecycle evidence", input, got)
			}
		})
	}
}

func TestExtractLifecycleEvidenceRejectsNilReader(t *testing.T) {
	got, err := ExtractLifecycleEvidence(nil)
	if err == nil {
		t.Fatalf("ExtractLifecycleEvidence(nil) returned nil error and evidence %#v", got)
	}
	if got != nil {
		t.Fatalf("evidence on error = %#v, want nil", got)
	}
	if !strings.Contains(err.Error(), "nil reader") {
		t.Fatalf("error = %q, want it to mention nil reader", err.Error())
	}
}

func TestLogOnlyEvidenceFromFailurePreservesFailureFieldsAndMarksLogSource(t *testing.T) {
	exit := 137
	failure := domain.LogOnlyFailure{
		ID:           "ab/c123def",
		Workdir:      "/runs/example/work/ab/c123def",
		Process:      "PIPE:ALIGN",
		Name:         "sample-01",
		Exit:         &exit,
		ErrorSummary: "killed by scheduler",
		ErrorBlock:   "ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
	}
	source := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log"}

	got := LogOnlyEvidenceFromFailure(failure, source)

	if got.ID != failure.ID || got.Workdir != failure.Workdir || got.Process != failure.Process || got.Name != failure.Name {
		t.Fatalf("identity fields = ID %q Workdir %q Process %q Name %q, want failure fields %#v", got.ID, got.Workdir, got.Process, got.Name, failure)
	}
	if got.ObservedStatus != domain.TaskStatusFailed {
		t.Fatalf("ObservedStatus = %q, want %q for failure-block evidence", got.ObservedStatus, domain.TaskStatusFailed)
	}
	if got.Exit == nil || *got.Exit != exit {
		if got.Exit == nil {
			t.Fatalf("Exit = nil, want %d", exit)
		}
		t.Fatalf("Exit = %d, want %d", *got.Exit, exit)
	}
	if got.ErrorSummary != failure.ErrorSummary || got.ErrorBlock != failure.ErrorBlock {
		t.Fatalf("error fields = Summary %q Block %q, want failure summary/block", got.ErrorSummary, got.ErrorBlock)
	}
	if got.Completeness != domain.LogOnlyEvidencePartial || got.CommandFilesAvailable {
		t.Fatalf("Completeness/CommandFilesAvailable = %q/%v, want partial/false for log-only failure evidence", got.Completeness, got.CommandFilesAvailable)
	}
	if len(got.Sources) != 1 {
		t.Fatalf("Sources = %#v, want one selected-log failure-block source", got.Sources)
	}
	logSource := got.Sources[0]
	if logSource.Kind != domain.LogOnlyEvidenceSourceLog || logSource.Path != source.Path {
		t.Fatalf("source = %#v, want log source at %q", logSource, source.Path)
	}
	if !strings.Contains(strings.ToLower(logSource.Detail), "failure") {
		t.Fatalf("source detail = %q, want it to identify failure-block provenance", logSource.Detail)
	}
}

func TestLogOnlyEvidenceFromFailureKeepsPartialFailureWhenFieldsAreMissing(t *testing.T) {
	failure := domain.LogOnlyFailure{
		Process:      "PIPE:QC",
		ErrorSummary: "task failed before workdir was reported",
		ErrorBlock:   "ERROR ~ Error executing process > 'PIPE:QC'",
	}

	got := LogOnlyEvidenceFromFailure(failure, domain.SourceFingerprint{Kind: domain.SourceKindLog})

	if got.ID != "" || got.Workdir != "" || got.Name != "" {
		t.Fatalf("ID/Workdir/Name = %q/%q/%q, want empty partial fields preserved", got.ID, got.Workdir, got.Name)
	}
	if got.Process != failure.Process || got.ErrorSummary != failure.ErrorSummary || got.ErrorBlock != failure.ErrorBlock {
		t.Fatalf("preserved fields = %#v, want process/error fields from %#v", got, failure)
	}
	if got.Exit != nil {
		t.Fatalf("Exit = %d, want nil partial failure exit", *got.Exit)
	}
	if got.ObservedStatus != domain.TaskStatusFailed || got.Completeness != domain.LogOnlyEvidencePartial || got.CommandFilesAvailable {
		t.Fatalf("status/completeness/command availability = %q/%q/%v, want failed/partial/false", got.ObservedStatus, got.Completeness, got.CommandFilesAvailable)
	}
	if len(got.Sources) != 1 || got.Sources[0].Kind != domain.LogOnlyEvidenceSourceLog {
		t.Fatalf("Sources = %#v, want one log source marker even when path is unavailable", got.Sources)
	}
}

func TestEnrichLogOnlyEvidenceFromWorkdirsAddsExitAndCommandError(t *testing.T) {
	workdir := t.TempDir()
	if err := os.WriteFile(filepath.Join(workdir, ".exitcode"), []byte("137\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(.exitcode) returned error: %v", err)
	}
	commandErr := strings.Join([]string{
		"setup complete",
		"fatal ERROR writing sample",
		"cleanup after failure",
	}, "\n")
	if err := os.WriteFile(filepath.Join(workdir, string(domain.CommandFileErr)), []byte(commandErr), 0o644); err != nil {
		t.Fatalf("WriteFile(.command.err) returned error: %v", err)
	}
	if err := os.WriteFile(filepath.Join(workdir, string(domain.CommandFileLog)), []byte("fallback log should not replace stderr"), 0o644); err != nil {
		t.Fatalf("WriteFile(.command.log) returned error: %v", err)
	}

	input := []domain.LogOnlyTaskEvidence{
		{
			ID:             "ab/c123def",
			Workdir:        workdir,
			Process:        "PIPE:ALIGN",
			ObservedStatus: domain.TaskStatusFailed,
			Sources: []domain.LogOnlyEvidenceSource{
				{Kind: domain.LogOnlyEvidenceSourceLog, Path: filepath.Join(workdir, "..", "..", ".nextflow.log"), Detail: "failure block"},
			},
			Completeness: domain.LogOnlyEvidencePartial,
		},
	}

	got, err := EnrichLogOnlyEvidenceFromWorkdirs(context.Background(), input, 256)
	if err != nil {
		t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs returned %d rows, want 1: %#v", len(got), got)
	}
	if input[0].Exit != nil || input[0].ErrorBlock != "" || input[0].CommandFilesAvailable {
		t.Fatalf("input evidence was mutated: %#v", input[0])
	}

	enriched := got[0]
	if enriched.Exit == nil || *enriched.Exit != 137 {
		if enriched.Exit == nil {
			t.Fatalf("Exit = nil, want 137 from .exitcode")
		}
		t.Fatalf("Exit = %d, want 137 from .exitcode", *enriched.Exit)
	}
	if enriched.ErrorBlock != commandErr {
		t.Fatalf("ErrorBlock = %q, want bounded .command.err content %q", enriched.ErrorBlock, commandErr)
	}
	if enriched.ErrorSummary != commandErr {
		t.Fatalf("ErrorSummary = %q, want deterministic summary from .command.err", enriched.ErrorSummary)
	}
	if !enriched.CommandFilesAvailable {
		t.Fatalf("CommandFilesAvailable = false, want true when referenced command evidence exists")
	}
	if enriched.Completeness != domain.LogOnlyEvidencePartial {
		t.Fatalf("Completeness = %q, want partial", enriched.Completeness)
	}

	exitSourcePath := filepath.Join(workdir, ".exitcode")
	errSourcePath := filepath.Join(workdir, string(domain.CommandFileErr))
	for _, want := range []string{exitSourcePath, errSourcePath} {
		found := false
		for _, source := range enriched.Sources {
			if source.Kind == domain.LogOnlyEvidenceSourceCommand && source.Path == want {
				found = true
				break
			}
		}
		if !found {
			t.Fatalf("Sources = %#v, want command-file source for %q", enriched.Sources, want)
		}
	}
}

func TestEnrichLogOnlyEvidenceFromWorkdirsUsesCommandLogFallbackWithByteBound(t *testing.T) {
	workdir := t.TempDir()
	commandLog := "ERROR long command log content that must be truncated before the end"
	if err := os.WriteFile(filepath.Join(workdir, string(domain.CommandFileLog)), []byte(commandLog), 0o644); err != nil {
		t.Fatalf("WriteFile(.command.log) returned error: %v", err)
	}

	got, err := EnrichLogOnlyEvidenceFromWorkdirs(context.Background(), []domain.LogOnlyTaskEvidence{{Workdir: workdir}}, 18)
	if err != nil {
		t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs returned %d rows, want 1", len(got))
	}
	if got[0].ErrorBlock != commandLog[:18] {
		t.Fatalf("ErrorBlock = %q, want first 18 bytes of .command.log %q", got[0].ErrorBlock, commandLog[:18])
	}
	if len(got[0].ErrorBlock) > 18 {
		t.Fatalf("ErrorBlock length = %d, want <= 18 bytes", len(got[0].ErrorBlock))
	}
	if got[0].ErrorSummary != commandLog[:18] {
		t.Fatalf("ErrorSummary = %q, want bounded summary %q", got[0].ErrorSummary, commandLog[:18])
	}
	if !got[0].CommandFilesAvailable {
		t.Fatalf("CommandFilesAvailable = false, want true for .command.log fallback")
	}

	foundLogSource := false
	for _, source := range got[0].Sources {
		if source.Kind == domain.LogOnlyEvidenceSourceCommand && source.Path == filepath.Join(workdir, string(domain.CommandFileLog)) && strings.Contains(source.Detail, "truncated") {
			foundLogSource = true
			break
		}
	}
	if !foundLogSource {
		t.Fatalf("Sources = %#v, want truncated command-file source for .command.log", got[0].Sources)
	}
}

func TestEnrichLogOnlyEvidenceFromWorkdirsSkipsMissingAndUnreferencedWorkdirs(t *testing.T) {
	runPath := t.TempDir()
	referencedMissing := filepath.Join(runPath, "work", "aa", "111111")
	unreferencedWorkdir := filepath.Join(runPath, "work", "bb", "222222")
	if err := os.MkdirAll(unreferencedWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(unreferencedWorkdir) returned error: %v", err)
	}
	if err := os.WriteFile(filepath.Join(unreferencedWorkdir, ".exitcode"), []byte("99\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(unreferenced .exitcode) returned error: %v", err)
	}
	if err := os.WriteFile(filepath.Join(unreferencedWorkdir, string(domain.CommandFileErr)), []byte("unreferenced error"), 0o644); err != nil {
		t.Fatalf("WriteFile(unreferenced .command.err) returned error: %v", err)
	}

	got, err := EnrichLogOnlyEvidenceFromWorkdirs(context.Background(), []domain.LogOnlyTaskEvidence{
		{ID: "aa/111111", Workdir: referencedMissing, Completeness: domain.LogOnlyEvidencePartial},
		{ID: "", Workdir: "", Process: "PIPE:NO_WORKDIR", Completeness: domain.LogOnlyEvidencePartial},
	}, 64)
	if err != nil {
		t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs returned error: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs returned %d rows, want 2: %#v", len(got), got)
	}
	for index, evidence := range got {
		if evidence.Exit != nil || evidence.ErrorSummary != "" || evidence.ErrorBlock != "" || evidence.CommandFilesAvailable {
			t.Fatalf("evidence[%d] = %#v, want missing/empty workdir evidence preserved without scanning unreferenced siblings", index, evidence)
		}
	}
}

func TestEnrichLogOnlyEvidenceFromWorkdirsRejectsInvalidInputs(t *testing.T) {
	workdir := t.TempDir()
	if err := os.WriteFile(filepath.Join(workdir, ".exitcode"), []byte("not-an-int\n"), 0o644); err != nil {
		t.Fatalf("WriteFile(.exitcode) returned error: %v", err)
	}

	t.Run("nil context", func(t *testing.T) {
		got, err := EnrichLogOnlyEvidenceFromWorkdirs(nil, []domain.LogOnlyTaskEvidence{{Workdir: workdir}}, 64)
		if err == nil {
			t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs(nil context) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on error = %#v, want nil", got)
		}
		if !strings.Contains(err.Error(), "nil context") {
			t.Fatalf("error = %q, want it to mention nil context", err.Error())
		}
	})

	t.Run("canceled context", func(t *testing.T) {
		ctx, cancel := context.WithCancel(context.Background())
		cancel()

		got, err := EnrichLogOnlyEvidenceFromWorkdirs(ctx, []domain.LogOnlyTaskEvidence{{Workdir: workdir}}, 64)
		if err == nil {
			t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs(canceled context) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on error = %#v, want nil", got)
		}
		if !errors.Is(err, context.Canceled) {
			t.Fatalf("error = %v, want context.Canceled", err)
		}
	})

	t.Run("non-positive byte bound", func(t *testing.T) {
		for _, maxBytes := range []int64{0, -1} {
			got, err := EnrichLogOnlyEvidenceFromWorkdirs(context.Background(), []domain.LogOnlyTaskEvidence{{Workdir: workdir}}, maxBytes)
			if err == nil {
				t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs(maxBytes=%d) returned nil error and evidence %#v", maxBytes, got)
			}
			if got != nil {
				t.Fatalf("evidence on maxBytes=%d error = %#v, want nil", maxBytes, got)
			}
			if !strings.Contains(err.Error(), "max bytes") {
				t.Fatalf("error = %q, want it to mention max bytes", err.Error())
			}
		}
	})

	t.Run("invalid exitcode", func(t *testing.T) {
		got, err := EnrichLogOnlyEvidenceFromWorkdirs(context.Background(), []domain.LogOnlyTaskEvidence{{Workdir: workdir}}, 64)
		if err == nil {
			t.Fatalf("EnrichLogOnlyEvidenceFromWorkdirs(invalid .exitcode) returned nil error and evidence %#v", got)
		}
		if got != nil {
			t.Fatalf("evidence on invalid .exitcode error = %#v, want nil", got)
		}
		for _, want := range []string{".exitcode", "invalid exit", "not-an-int"} {
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
			}
		}
	})
}

func TestExtractFailureBlocksParsesCommonNextflowErrorBlock(t *testing.T) {
	input := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:01:00.000 [Task monitor] ERROR nextflow.processor.TaskProcessor - Error executing process > 'PIPE:ALIGN (sample-01)'",
		"",
		"Caused by:",
		"  Process `PIPE:ALIGN (sample-01)` terminated with an error exit status (137)",
		"",
		"Command exit status:",
		"  137",
		"",
		"Command error:",
		"  killed by scheduler",
		"",
		"Work dir:",
		"  /runs/example/work/ab/c123def",
		"",
		"Tip: you can replicate the issue by changing to the process work dir and entering `bash .command.run`",
		"Apr-28 12:01:01.000 [main] INFO nextflow.Session - Execution stopped",
		"",
	}, "\n")

	got, err := ExtractFailureBlocks(strings.NewReader(input))
	if err != nil {
		t.Fatalf("ExtractFailureBlocks returned error: %v", err)
	}
	if len(got) != 1 {
		t.Fatalf("ExtractFailureBlocks returned %d blocks, want 1: %#v", len(got), got)
	}

	block := got[0]
	if block.Process != "PIPE:ALIGN" {
		t.Fatalf("Process = %q, want %q", block.Process, "PIPE:ALIGN")
	}
	if block.Name != "sample-01" {
		t.Fatalf("Name = %q, want %q", block.Name, "sample-01")
	}
	if block.Workdir != "/runs/example/work/ab/c123def" {
		t.Fatalf("Workdir = %q, want %q", block.Workdir, "/runs/example/work/ab/c123def")
	}
	if block.Exit == nil || *block.Exit != 137 {
		if block.Exit == nil {
			t.Fatalf("Exit = nil, want 137")
		}
		t.Fatalf("Exit = %d, want 137", *block.Exit)
	}

	for _, want := range []string{"Error executing process", "Command error:", "Work dir:"} {
		if !strings.Contains(block.Block, want) {
			t.Fatalf("Block = %q, want it to contain %q", block.Block, want)
		}
	}
	for _, notWant := range []string{"Session start", "Execution stopped"} {
		if strings.Contains(block.Block, notWant) {
			t.Fatalf("Block = %q, want it not to contain surrounding log line %q", block.Block, notWant)
		}
	}
}

func TestExtractFailureBlocksParsesMultipleBlocksInOrderWithPartialEvidence(t *testing.T) {
	input := strings.Join([]string{
		"ERROR ~ Error executing process > 'PIPE:QC'",
		"",
		"Caused by:",
		"  Process `PIPE:QC` terminated with an error exit status (1)",
		"",
		"Work dir:",
		"  work/de/f456",
		"Apr-28 12:02:00.000 [main] DEBUG nextflow.Session - cleanup after first failure",
		"Apr-28 12:03:00.000 [Task monitor] ERROR nextflow.processor.TaskProcessor - Process `PIPE:CALL (tumor-02)` terminated with an error exit status (2)",
		"Apr-28 12:03:01.000 [main] INFO nextflow.Session - Execution stopped",
		"",
	}, "\n")

	got, err := ExtractFailureBlocks(strings.NewReader(input))
	if err != nil {
		t.Fatalf("ExtractFailureBlocks returned error: %v", err)
	}
	if len(got) != 2 {
		t.Fatalf("ExtractFailureBlocks returned %d blocks, want 2: %#v", len(got), got)
	}

	if got[0].Process != "PIPE:QC" || got[0].Name != "" || got[0].Workdir != "work/de/f456" {
		t.Fatalf("first block fields = Process %q, Name %q, Workdir %q; want PIPE:QC, empty name, work/de/f456", got[0].Process, got[0].Name, got[0].Workdir)
	}
	if got[0].Exit == nil || *got[0].Exit != 1 {
		if got[0].Exit == nil {
			t.Fatalf("first block Exit = nil, want 1")
		}
		t.Fatalf("first block Exit = %d, want 1", *got[0].Exit)
	}

	if got[1].Process != "PIPE:CALL" || got[1].Name != "tumor-02" || got[1].Workdir != "" {
		t.Fatalf("second block fields = Process %q, Name %q, Workdir %q; want PIPE:CALL, tumor-02, empty workdir", got[1].Process, got[1].Name, got[1].Workdir)
	}
	if got[1].Exit == nil || *got[1].Exit != 2 {
		if got[1].Exit == nil {
			t.Fatalf("second block Exit = nil, want 2")
		}
		t.Fatalf("second block Exit = %d, want 2", *got[1].Exit)
	}
	if strings.Contains(got[0].Block, "cleanup after first failure") {
		t.Fatalf("first block = %q, want next timestamped log line excluded", got[0].Block)
	}
}

func TestExtractFailureBlocksReturnsEmptySliceWhenNoFailureEvidenceAppears(t *testing.T) {
	input := strings.Join([]string{
		"Apr-28 12:00:00.000 [main] INFO nextflow.Session - Session start",
		"Apr-28 12:01:00.000 [main] WARN nextflow.Session - A warning without task failure evidence",
		"Apr-28 12:02:00.000 [main] ERROR nextflow.Session - Pipeline aborted before any process failure was reported",
		"",
	}, "\n")

	got, err := ExtractFailureBlocks(strings.NewReader(input))
	if err != nil {
		t.Fatalf("ExtractFailureBlocks returned error: %v", err)
	}
	if len(got) != 0 {
		t.Fatalf("ExtractFailureBlocks returned %#v, want no blocks", got)
	}
}

func TestExtractFailureBlocksRejectsNilReader(t *testing.T) {
	got, err := ExtractFailureBlocks(nil)
	if err == nil {
		t.Fatalf("ExtractFailureBlocks(nil) returned nil error and blocks %#v", got)
	}
	if !strings.Contains(err.Error(), "nil reader") {
		t.Fatalf("error = %q, want it to mention nil reader", err.Error())
	}
}

func TestExtractWorkdirEvidenceParsesCommonForms(t *testing.T) {
	tests := []struct {
		name        string
		block       string
		wantID      string
		wantWorkdir string
	}{
		{
			name: "work dir on following line with absolute path",
			block: strings.Join([]string{
				"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
				"Work dir:",
				"  /runs/example/work/ab/c123def",
			}, "\n"),
			wantID:      "ab/c123def",
			wantWorkdir: "/runs/example/work/ab/c123def",
		},
		{
			name:        "working directory on same line with relative work path",
			block:       "Working directory: `work/DE/F456`;",
			wantID:      "de/f456",
			wantWorkdir: "work/DE/F456",
		},
		{
			name:        "hash-only work dir value",
			block:       "Work dir: [AB/C123DEF]",
			wantID:      "ab/c123def",
			wantWorkdir: "",
		},
		{
			name:        "inline nested command file below work path",
			block:       "See /runs/example/work/12/abcdef/.command.err for command stderr.",
			wantID:      "12/abcdef",
			wantWorkdir: "/runs/example/work/12/abcdef",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotID, gotWorkdir, err := ExtractWorkdirEvidence(tt.block)
			if err != nil {
				t.Fatalf("ExtractWorkdirEvidence(%q) returned error: %v", tt.block, err)
			}
			if gotID != tt.wantID || gotWorkdir != tt.wantWorkdir {
				t.Fatalf("ExtractWorkdirEvidence(%q) = (%q, %q), want (%q, %q)", tt.block, gotID, gotWorkdir, tt.wantID, tt.wantWorkdir)
			}
		})
	}
}

func TestExtractWorkdirEvidenceReturnsEmptyWhenMissing(t *testing.T) {
	blocks := []string{
		"",
		"  \t\n  ",
		strings.Join([]string{
			"ERROR ~ Error executing process > 'PIPE:QC'",
			"Command error:",
			"  killed by scheduler before workdir was reported",
			"Tip: you can replicate the issue by changing to the process work dir and entering `bash .command.run`",
		}, "\n"),
	}

	for _, block := range blocks {
		t.Run(block, func(t *testing.T) {
			gotID, gotWorkdir, err := ExtractWorkdirEvidence(block)
			if err != nil {
				t.Fatalf("ExtractWorkdirEvidence(%q) returned error: %v", block, err)
			}
			if gotID != "" || gotWorkdir != "" {
				t.Fatalf("ExtractWorkdirEvidence(%q) = (%q, %q), want empty id and workdir", block, gotID, gotWorkdir)
			}
		})
	}
}

func TestExtractWorkdirEvidenceRejectsInvalidDirectWorkdirLabel(t *testing.T) {
	tests := []string{
		"Work dir: not-a-workdir",
		strings.Join([]string{
			"Work dir:",
			"  /runs/example/work/zz/c123def",
		}, "\n"),
	}

	for _, block := range tests {
		t.Run(block, func(t *testing.T) {
			gotID, gotWorkdir, err := ExtractWorkdirEvidence(block)
			if err == nil {
				t.Fatalf("ExtractWorkdirEvidence(%q) returned nil error and values (%q, %q)", block, gotID, gotWorkdir)
			}
			if gotID != "" || gotWorkdir != "" {
				t.Fatalf("values on error = (%q, %q), want empty id and workdir", gotID, gotWorkdir)
			}
			for _, want := range []string{"extract workdir evidence", "invalid workdir"} {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestExtractExitEvidenceParsesCommonForms(t *testing.T) {
	tests := []struct {
		name  string
		block string
		want  int
	}{
		{
			name:  "command exit status on same line",
			block: "Command exit status: 1",
			want:  1,
		},
		{
			name:  "lowercase exit status on same line",
			block: "exit status: 2",
			want:  2,
		},
		{
			name:  "uppercase exit status on same line",
			block: "Exit status: 137",
			want:  137,
		},
		{
			name: "command exit status on following line",
			block: strings.Join([]string{
				"Command exit status:",
				"  0",
				"Command error:",
				"  no error text captured",
			}, "\n"),
			want: 0,
		},
		{
			name:  "nextflow error exit status parenthetical",
			block: "Process `PIPE:QC` terminated with an error exit status (143)",
			want:  143,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ExtractExitEvidence(tt.block)
			if err != nil {
				t.Fatalf("ExtractExitEvidence(%q) returned error: %v", tt.block, err)
			}
			if got == nil {
				t.Fatalf("ExtractExitEvidence(%q) = nil, want %d", tt.block, tt.want)
			}
			if *got != tt.want {
				t.Fatalf("ExtractExitEvidence(%q) = %d, want %d", tt.block, *got, tt.want)
			}
		})
	}
}

func TestExtractExitEvidenceReturnsNilWhenMissing(t *testing.T) {
	block := strings.Join([]string{
		"ERROR ~ Error executing process > 'PIPE:QC'",
		"Command error:",
		"  killed by scheduler before an exit line was reported",
		"Work dir:",
		"  /runs/example/work/ab/c123def",
	}, "\n")

	got, err := ExtractExitEvidence(block)
	if err != nil {
		t.Fatalf("ExtractExitEvidence returned error: %v", err)
	}
	if got != nil {
		t.Fatalf("ExtractExitEvidence = %d, want nil", *got)
	}
}

func TestExtractExitEvidenceRejectsInvalidDirectlyLabelledExit(t *testing.T) {
	tests := []string{
		"Command exit status: killed by scheduler",
		"Exit status: one",
		strings.Join([]string{
			"Command exit status:",
			"  failed",
		}, "\n"),
	}

	for _, block := range tests {
		t.Run(block, func(t *testing.T) {
			got, err := ExtractExitEvidence(block)
			if err == nil {
				t.Fatalf("ExtractExitEvidence(%q) returned nil error and value %#v", block, got)
			}
			if got != nil {
				t.Fatalf("value on error = %d, want nil", *got)
			}
			for _, want := range []string{"extract exit evidence", "invalid exit"} {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestSummarizeErrorBlockReturnsEmptyForBlankBlock(t *testing.T) {
	for _, block := range []string{"", "   \n\t\n", "\r\n\r\n"} {
		if got := SummarizeErrorBlock(block, 80); got != "" {
			t.Fatalf("SummarizeErrorBlock(%q, 80) = %q, want empty summary", block, got)
		}
	}
}

func TestSummarizeErrorBlockKeepsShortCommandErrorText(t *testing.T) {
	block := strings.Join([]string{
		"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
		"",
		"Command error:",
		"  bwa: failed to open input.bam",
		"  exit status 2",
		"",
		"Work dir:",
		"  /runs/example/work/ab/c123def",
		"",
		"Tip: you can replicate the issue by changing to the process work dir and entering `bash .command.run`",
	}, "\n")

	want := strings.Join([]string{
		"bwa: failed to open input.bam",
		"exit status 2",
	}, "\n")
	if got := SummarizeErrorBlock(block, 200); got != want {
		t.Fatalf("SummarizeErrorBlock returned %q, want %q", got, want)
	}
}

func TestSummarizeErrorBlockTruncatesLongSummaryToMaxBytes(t *testing.T) {
	block := strings.Join([]string{
		"Command error:",
		"  0123456789",
		"  abcdefghij",
		"  klmnopqrst",
		"Work dir:",
		"  /runs/example/work/ab/c123def",
	}, "\n")

	const maxBytes = 15
	got := SummarizeErrorBlock(block, maxBytes)
	want := "0123456789\nabcd"
	if got != want {
		t.Fatalf("SummarizeErrorBlock returned %q, want %q", got, want)
	}
	if len(got) > maxBytes {
		t.Fatalf("summary length = %d bytes, want at most %d", len(got), maxBytes)
	}
}

func TestNormalizeFailureBlockUsesParsedWorkdirAndExitBeforeTextFallback(t *testing.T) {
	runPath := t.TempDir()
	derivedWorkdir := filepath.Join(runPath, "work", "ab", "c123def")
	if err := os.MkdirAll(derivedWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", derivedWorkdir, err)
	}

	parsedExit := 137
	block := FailureBlock{
		Process: "PIPE:ALIGN",
		Name:    "sample-01",
		Workdir: "AB/C123DEF",
		Exit:    &parsedExit,
		Block: strings.Join([]string{
			"ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
			"Command exit status: 1",
			"Command error:",
			"  parsed fields should win over text fallback",
			"Work dir:",
			"  /other/run/work/de/f456",
		}, "\n"),
	}

	got, err := NormalizeFailureBlock(domain.RunDir{Path: runPath}, block)
	if err != nil {
		t.Fatalf("NormalizeFailureBlock returned error: %v", err)
	}

	if got.ID != "ab/c123def" {
		t.Fatalf("ID = %q, want %q", got.ID, "ab/c123def")
	}
	if got.Workdir != derivedWorkdir {
		t.Fatalf("Workdir = %q, want derived existing path %q", got.Workdir, derivedWorkdir)
	}
	if got.Process != block.Process || got.Name != block.Name {
		t.Fatalf("Process/Name = %q/%q, want %q/%q", got.Process, got.Name, block.Process, block.Name)
	}
	if got.Exit == nil || *got.Exit != parsedExit {
		if got.Exit == nil {
			t.Fatalf("Exit = nil, want %d", parsedExit)
		}
		t.Fatalf("Exit = %d, want parsed field %d", *got.Exit, parsedExit)
	}
	if got.ErrorSummary != "parsed fields should win over text fallback" {
		t.Fatalf("ErrorSummary = %q, want command error summary", got.ErrorSummary)
	}
	if got.ErrorBlock != block.Block {
		t.Fatalf("ErrorBlock = %q, want original block %q", got.ErrorBlock, block.Block)
	}
}

func TestNormalizeFailureBlockExtractsMissingWorkdirAndExitFromBlockText(t *testing.T) {
	runPath := t.TempDir()
	derivedWorkdir := filepath.Join(runPath, "work", "bb", "222222")
	if err := os.MkdirAll(derivedWorkdir, 0o755); err != nil {
		t.Fatalf("MkdirAll(%q) returned error: %v", derivedWorkdir, err)
	}

	block := FailureBlock{
		Process: "PIPE:CALL",
		Name:    "tumor-02",
		Block: strings.Join([]string{
			"ERROR ~ Error executing process > 'PIPE:CALL (tumor-02)'",
			"Command exit status:",
			"  2",
			"Command error:",
			"  No such file or directory",
			"Work dir:",
			"  BB/222222",
		}, "\n"),
	}

	got, err := NormalizeFailureBlock(domain.RunDir{Path: runPath}, block)
	if err != nil {
		t.Fatalf("NormalizeFailureBlock returned error: %v", err)
	}

	if got.ID != "bb/222222" {
		t.Fatalf("ID = %q, want %q", got.ID, "bb/222222")
	}
	if got.Workdir != derivedWorkdir {
		t.Fatalf("Workdir = %q, want derived existing path %q", got.Workdir, derivedWorkdir)
	}
	if got.Exit == nil || *got.Exit != 2 {
		if got.Exit == nil {
			t.Fatalf("Exit = nil, want 2")
		}
		t.Fatalf("Exit = %d, want 2", *got.Exit)
	}
	if got.ErrorSummary != "No such file or directory" {
		t.Fatalf("ErrorSummary = %q, want %q", got.ErrorSummary, "No such file or directory")
	}
	if got.Process != block.Process || got.Name != block.Name {
		t.Fatalf("Process/Name = %q/%q, want %q/%q", got.Process, got.Name, block.Process, block.Name)
	}
}

func TestNormalizeFailureBlockAllowsPartialEvidenceWhenWorkdirAndExitAreMissing(t *testing.T) {
	block := FailureBlock{
		Process: "PIPE:QC",
		Block: strings.Join([]string{
			"ERROR ~ Error executing process > 'PIPE:QC'",
			"Command error:",
			"  scheduler killed task before metadata was written",
		}, "\n"),
	}

	got, err := NormalizeFailureBlock(domain.RunDir{Path: t.TempDir()}, block)
	if err != nil {
		t.Fatalf("NormalizeFailureBlock returned error: %v", err)
	}

	if got.ID != "" || got.Workdir != "" {
		t.Fatalf("ID/Workdir = %q/%q, want empty partial evidence", got.ID, got.Workdir)
	}
	if got.Exit != nil {
		t.Fatalf("Exit = %d, want nil partial evidence", *got.Exit)
	}
	if got.Process != "PIPE:QC" || got.Name != "" {
		t.Fatalf("Process/Name = %q/%q, want PIPE:QC/empty", got.Process, got.Name)
	}
	if got.ErrorSummary != "scheduler killed task before metadata was written" {
		t.Fatalf("ErrorSummary = %q, want deterministic command error summary", got.ErrorSummary)
	}
}

func TestNormalizeFailureBlockReturnsEvidenceExtractionErrors(t *testing.T) {
	tests := []struct {
		name      string
		block     FailureBlock
		wantParts []string
	}{
		{
			name: "invalid workdir in parsed field",
			block: FailureBlock{
				Workdir: "not-a-workdir",
				Block:   "Work dir: /runs/example/work/ab/c123def",
			},
			wantParts: []string{"normalize failure block", "parsed workdir", "not-a-workdir"},
		},
		{
			name:      "invalid workdir in block text",
			block:     FailureBlock{Block: "Work dir: not-a-workdir"},
			wantParts: []string{"normalize failure block", "invalid workdir"},
		},
		{
			name:      "invalid exit in block text",
			block:     FailureBlock{Block: "Command exit status: failed"},
			wantParts: []string{"normalize failure block", "invalid exit"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := NormalizeFailureBlock(domain.RunDir{Path: t.TempDir()}, tt.block)
			if err == nil {
				t.Fatalf("NormalizeFailureBlock returned nil error and evidence %#v", got)
			}
			for _, want := range tt.wantParts {
				if !strings.Contains(err.Error(), want) {
					t.Fatalf("error = %q, want it to mention %q", err.Error(), want)
				}
			}
		})
	}
}

func TestBuildLogOnlyEvidenceStatusSummarizesObservedEvidenceAsIncomplete(t *testing.T) {
	exit := 2
	runDir := domain.RunDir{Path: "/runs/log-only"}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/log-only/.nextflow.log", Size: 4096}
	artifacts := domain.ArtifactSet{
		RunDir:           runDir,
		Mode:             domain.IndexModeLogOnly,
		Log:              &logSource,
		SearchedPatterns: []string{"trace*.txt", ".nextflow*.log"},
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticInfo, Code: "selected_log", Message: "selected newest Nextflow log"},
		},
	}
	evidence := []domain.LogOnlyTaskEvidence{
		{
			ID:             "ab/c123def",
			Process:        "PIPE:ALIGN",
			Name:           "sample-01",
			ObservedStatus: domain.TaskStatusCompleted,
			Sources:        []domain.LogOnlyEvidenceSource{{Kind: domain.LogOnlyEvidenceSourceLog, Path: logSource.Path}},
			Completeness:   domain.LogOnlyEvidencePartial,
		},
		{
			ID:                    "bb/222222",
			Workdir:               "/runs/log-only/work/bb/222222",
			Process:               "PIPE:CALL",
			Name:                  "tumor-02",
			ObservedStatus:        domain.TaskStatusFailed,
			Exit:                  &exit,
			ErrorSummary:          "No such file or directory",
			ErrorBlock:            "ERROR ~ Error executing process > 'PIPE:CALL (tumor-02)'",
			Sources:               []domain.LogOnlyEvidenceSource{{Kind: domain.LogOnlyEvidenceSourceLog, Path: logSource.Path}},
			Completeness:          domain.LogOnlyEvidencePartial,
			CommandFilesAvailable: true,
		},
		{
			ID:             "cc/333333",
			Process:        "PIPE:STOP",
			ObservedStatus: domain.TaskStatusAborted,
			ErrorSummary:   "workflow aborted",
			Completeness:   domain.LogOnlyEvidencePartial,
		},
		{
			ID:             "dd/444444",
			Process:        "PIPE:CACHE",
			ObservedStatus: domain.TaskStatusCached,
			Completeness:   domain.LogOnlyEvidencePartial,
		},
	}

	got, err := BuildLogOnlyEvidenceStatus(runDir, artifacts, evidence)
	if err != nil {
		t.Fatalf("BuildLogOnlyEvidenceStatus returned error: %v", err)
	}

	if got.RunDir != runDir {
		t.Fatalf("RunDir = %#v, want %#v", got.RunDir, runDir)
	}
	if got.Mode != domain.IndexModeLogOnly {
		t.Fatalf("Mode = %q, want %q", got.Mode, domain.IndexModeLogOnly)
	}
	if got.Freshness != domain.IndexFreshnessUnsupported {
		t.Fatalf("Freshness = %q, want %q for degraded log-only evidence", got.Freshness, domain.IndexFreshnessUnsupported)
	}
	if got.IndexPath != "" || got.BuiltAt != nil {
		t.Fatalf("IndexPath/BuiltAt = %q/%#v, want no trace-backed index metadata", got.IndexPath, got.BuiltAt)
	}
	if got.Sources.Log == nil || *got.Sources.Log != logSource {
		t.Fatalf("Sources.Log = %#v, want %#v", got.Sources.Log, logSource)
	}

	assertStatusCount(t, got.Counts, domain.TaskStatusAborted, 1)
	assertStatusCount(t, got.Counts, domain.TaskStatusCached, 1)
	assertStatusCount(t, got.Counts, domain.TaskStatusCompleted, 1)
	assertStatusCount(t, got.Counts, domain.TaskStatusFailed, 1)
	if len(got.Counts) != 4 {
		t.Fatalf("Counts = %#v, want exactly four observed status counts", got.Counts)
	}
	if got.FailedCount != 2 {
		t.Fatalf("FailedCount = %d, want 2 observed failed/aborted evidence rows", got.FailedCount)
	}
	if len(got.FailedPreview) != 0 {
		t.Fatalf("FailedPreview = %#v, want empty because log-only evidence is not trace-backed task rows", got.FailedPreview)
	}
	if len(got.LogOnlyEvidence) != len(evidence) {
		t.Fatalf("LogOnlyEvidence length = %d, want %d", len(got.LogOnlyEvidence), len(evidence))
	}
	if got.LogOnlyEvidence[1].ID != evidence[1].ID || got.LogOnlyEvidence[1].Exit == nil || *got.LogOnlyEvidence[1].Exit != exit || !got.LogOnlyEvidence[1].CommandFilesAvailable {
		t.Fatalf("LogOnlyEvidence[1] = %#v, want failed evidence fields preserved", got.LogOnlyEvidence[1])
	}
	got.LogOnlyEvidence[0].Process = "mutated output"
	if evidence[0].Process != "PIPE:ALIGN" {
		t.Fatalf("input evidence was aliased and mutated: %#v", evidence[0])
	}

	assertDiagnosticWithCode(t, got.Diagnostics, "selected_log")
	degraded := assertDiagnosticWithCode(t, got.Diagnostics, "log_only_degraded")
	if degraded.Severity != domain.DiagnosticWarning {
		t.Fatalf("log_only_degraded severity = %q, want %q", degraded.Severity, domain.DiagnosticWarning)
	}
	for _, want := range []string{"observed", "incomplete", "trace", logSource.Path, "trace*.txt"} {
		if !strings.Contains(degraded.Message+"\n"+degraded.Detail, want) {
			t.Fatalf("log_only_degraded diagnostic = %#v, want it to mention %q", degraded, want)
		}
	}

	recommendation := assertDiagnosticWithCode(t, got.Diagnostics, "nextflow_with_trace_recommended")
	if recommendation.Severity != domain.DiagnosticInfo {
		t.Fatalf("nextflow_with_trace_recommended severity = %q, want %q", recommendation.Severity, domain.DiagnosticInfo)
	}
	if !strings.Contains(recommendation.Detail, "-with-trace") {
		t.Fatalf("nextflow_with_trace_recommended detail = %q, want it to mention -with-trace", recommendation.Detail)
	}
}

func TestBuildLogOnlyEvidenceStatusExplainsNoObservedEvidence(t *testing.T) {
	runDir := domain.RunDir{Path: "/runs/no-trace"}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/no-trace/.nextflow.log", Size: 2048}
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeLogOnly, Log: &logSource}

	got, err := BuildLogOnlyEvidenceStatus(runDir, artifacts, nil)
	if err != nil {
		t.Fatalf("BuildLogOnlyEvidenceStatus returned error for empty evidence: %v", err)
	}

	if got.Mode != domain.IndexModeLogOnly {
		t.Fatalf("Mode = %q, want %q", got.Mode, domain.IndexModeLogOnly)
	}
	if len(got.Counts) != 0 {
		t.Fatalf("Counts = %#v, want no observed counts when no log-only evidence was parsed", got.Counts)
	}
	if got.FailedCount != 0 {
		t.Fatalf("FailedCount = %d, want 0 observed failures", got.FailedCount)
	}
	if len(got.LogOnlyEvidence) != 0 {
		t.Fatalf("LogOnlyEvidence = %#v, want empty", got.LogOnlyEvidence)
	}
	if len(got.FailedPreview) != 0 || len(got.LogOnlyFailures) != 0 {
		t.Fatalf("FailedPreview/LogOnlyFailures = %#v/%#v, want neither for no-evidence log-only status", got.FailedPreview, got.LogOnlyFailures)
	}

	missingEvidence := assertDiagnosticWithCode(t, got.Diagnostics, "log_only_no_parseable_evidence")
	if missingEvidence.Severity != domain.DiagnosticError {
		t.Fatalf("log_only_no_parseable_evidence severity = %q, want %q", missingEvidence.Severity, domain.DiagnosticError)
	}
	for _, want := range []string{"No parseable", "task", logSource.Path, "complete task"} {
		if !strings.Contains(missingEvidence.Message+"\n"+missingEvidence.Detail, want) {
			t.Fatalf("log_only_no_parseable_evidence diagnostic = %#v, want it to mention %q", missingEvidence, want)
		}
	}

	recommendation := assertDiagnosticWithCode(t, got.Diagnostics, "nextflow_with_trace_recommended")
	if recommendation.Severity != domain.DiagnosticInfo {
		t.Fatalf("nextflow_with_trace_recommended severity = %q, want %q", recommendation.Severity, domain.DiagnosticInfo)
	}
	if !strings.Contains(recommendation.Detail, "-with-trace") {
		t.Fatalf("nextflow_with_trace_recommended detail = %q, want it to mention -with-trace", recommendation.Detail)
	}
}

func TestBuildLogOnlyStatusReturnsDegradedSummaryWithFailureEvidence(t *testing.T) {
	exit := 137
	runDir := domain.RunDir{Path: "/runs/example"}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/example/.nextflow.log", Size: 4096}
	artifacts := domain.ArtifactSet{
		RunDir:           runDir,
		Mode:             domain.IndexModeLogOnly,
		Log:              &logSource,
		SearchedPatterns: []string{".nextflow.log"},
		Diagnostics: []domain.Diagnostic{
			{Severity: domain.DiagnosticInfo, Code: "selected_log", Message: "selected newest Nextflow log"},
		},
	}
	failures := []domain.LogOnlyFailure{
		{
			ID:           "ab/c123def",
			Workdir:      "/runs/example/work/ab/c123def",
			Process:      "PIPE:ALIGN",
			Name:         "sample-01",
			Exit:         &exit,
			ErrorSummary: "killed by scheduler",
			ErrorBlock:   "ERROR ~ Error executing process > 'PIPE:ALIGN (sample-01)'",
		},
	}

	got, err := BuildLogOnlyStatus(runDir, artifacts, failures)
	if err != nil {
		t.Fatalf("BuildLogOnlyStatus returned error: %v", err)
	}

	if got.RunDir != runDir {
		t.Fatalf("RunDir = %#v, want %#v", got.RunDir, runDir)
	}
	if got.Mode != domain.IndexModeLogOnly {
		t.Fatalf("Mode = %q, want %q", got.Mode, domain.IndexModeLogOnly)
	}
	if got.Freshness != domain.IndexFreshnessUnsupported {
		t.Fatalf("Freshness = %q, want %q for non-trace-backed log-only status", got.Freshness, domain.IndexFreshnessUnsupported)
	}
	if got.IndexPath != "" {
		t.Fatalf("IndexPath = %q, want empty because log-only status has no complete index", got.IndexPath)
	}
	if got.BuiltAt != nil {
		t.Fatalf("BuiltAt = %#v, want nil because no complete task table was built", got.BuiltAt)
	}
	if got.Sources.Log == nil || *got.Sources.Log != logSource {
		t.Fatalf("Sources.Log = %#v, want %#v", got.Sources.Log, logSource)
	}
	if len(got.Counts) != 0 {
		t.Fatalf("Counts = %#v, want no status counts because log-only evidence is not a complete task table", got.Counts)
	}
	if len(got.FailedPreview) != 0 {
		t.Fatalf("FailedPreview = %#v, want empty because log-only evidence must not be presented as trace-backed task rows", got.FailedPreview)
	}
	if got.FailedCount != len(failures) {
		t.Fatalf("FailedCount = %d, want %d parseable log-only failures", got.FailedCount, len(failures))
	}
	if len(got.LogOnlyFailures) != len(failures) {
		t.Fatalf("LogOnlyFailures length = %d, want %d", len(got.LogOnlyFailures), len(failures))
	}
	if got.LogOnlyFailures[0].Process != failures[0].Process || got.LogOnlyFailures[0].Name != failures[0].Name || got.LogOnlyFailures[0].Workdir != failures[0].Workdir {
		t.Fatalf("LogOnlyFailures[0] = %#v, want evidence %#v", got.LogOnlyFailures[0], failures[0])
	}
	if got.LogOnlyFailures[0].Exit == nil || *got.LogOnlyFailures[0].Exit != exit {
		t.Fatalf("LogOnlyFailures[0].Exit = %#v, want %d", got.LogOnlyFailures[0].Exit, exit)
	}

	assertDiagnosticWithCode(t, got.Diagnostics, "selected_log")
	degraded := assertDiagnosticWithCode(t, got.Diagnostics, "log_only_degraded")
	if degraded.Severity != domain.DiagnosticWarning {
		t.Fatalf("log_only_degraded severity = %q, want %q", degraded.Severity, domain.DiagnosticWarning)
	}
	for _, want := range []string{"log-only", "complete task", "resource", "status"} {
		if !strings.Contains(degraded.Message+"\n"+degraded.Detail, want) {
			t.Fatalf("log_only_degraded diagnostic = %#v, want it to mention %q", degraded, want)
		}
	}
}

func TestBuildLogOnlyStatusExplainsNoParseableFailures(t *testing.T) {
	runDir := domain.RunDir{Path: "/runs/no-trace"}
	logSource := domain.SourceFingerprint{Kind: domain.SourceKindLog, Path: "/runs/no-trace/.nextflow.log", Size: 2048}
	artifacts := domain.ArtifactSet{RunDir: runDir, Mode: domain.IndexModeLogOnly, Log: &logSource}

	got, err := BuildLogOnlyStatus(runDir, artifacts, nil)
	if err != nil {
		t.Fatalf("BuildLogOnlyStatus returned error for empty failure evidence: %v", err)
	}

	if got.Mode != domain.IndexModeLogOnly {
		t.Fatalf("Mode = %q, want %q", got.Mode, domain.IndexModeLogOnly)
	}
	if got.FailedCount != 0 {
		t.Fatalf("FailedCount = %d, want 0 parseable log-only failures", got.FailedCount)
	}
	if len(got.LogOnlyFailures) != 0 {
		t.Fatalf("LogOnlyFailures = %#v, want empty", got.LogOnlyFailures)
	}
	if len(got.Counts) != 0 || len(got.FailedPreview) != 0 {
		t.Fatalf("Counts/FailedPreview = %#v/%#v, want neither for log-only no-evidence status", got.Counts, got.FailedPreview)
	}

	missingEvidence := assertDiagnosticWithCode(t, got.Diagnostics, "log_only_no_parseable_failures")
	if missingEvidence.Severity != domain.DiagnosticError {
		t.Fatalf("log_only_no_parseable_failures severity = %q, want %q", missingEvidence.Severity, domain.DiagnosticError)
	}
	for _, want := range []string{"No parseable", logSource.Path, "complete task"} {
		if !strings.Contains(missingEvidence.Message+"\n"+missingEvidence.Detail, want) {
			t.Fatalf("log_only_no_parseable_failures diagnostic = %#v, want it to mention %q", missingEvidence, want)
		}
	}

	recommendation := assertDiagnosticWithCode(t, got.Diagnostics, "nextflow_with_trace_recommended")
	if recommendation.Severity != domain.DiagnosticInfo {
		t.Fatalf("nextflow_with_trace_recommended severity = %q, want %q", recommendation.Severity, domain.DiagnosticInfo)
	}
	if !strings.Contains(recommendation.Detail, "-with-trace") {
		t.Fatalf("nextflow_with_trace_recommended detail = %q, want it to mention -with-trace", recommendation.Detail)
	}
}

func assertStatusCount(t *testing.T, counts []domain.StatusCount, status domain.TaskStatus, want int) {
	t.Helper()
	for _, count := range counts {
		if count.Status == status {
			if count.Count != want {
				t.Fatalf("count for status %q = %d, want %d in %#v", status, count.Count, want, counts)
			}
			return
		}
	}
	t.Fatalf("count for status %q not found in %#v", status, counts)
}

func assertDiagnosticWithCode(t *testing.T, diagnostics []domain.Diagnostic, code string) domain.Diagnostic {
	t.Helper()
	for _, diagnostic := range diagnostics {
		if diagnostic.Code == code {
			return diagnostic
		}
	}
	t.Fatalf("diagnostic with code %q not found in %#v", code, diagnostics)
	return domain.Diagnostic{}
}
