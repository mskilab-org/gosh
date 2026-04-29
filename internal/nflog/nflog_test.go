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
