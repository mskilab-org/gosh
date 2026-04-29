package docs

import (
	"os"
	"strings"
	"testing"
)

func readGoV1TriageDoc(t *testing.T) string {
	t.Helper()

	content, err := os.ReadFile("go-v1-triage.md")
	if err != nil {
		t.Fatalf("read go-v1-triage.md: %v", err)
	}
	return string(content)
}

func TestGoV1TriageDocsCoverImplementedScope(t *testing.T) {
	doc := readGoV1TriageDoc(t)

	if strings.Contains(strings.ToLower(doc), "todo") {
		t.Fatalf("go-v1-triage.md still contains TODO text")
	}

	required := []string{
		"# Go v1 fast triage CLI",
		"read-only",
		"does not execute Nextflow",
		"does not invoke Python",
		"`status`",
		"`tasks`",
		"`inspect`",
		"`index`",
		"`--run-dir`",
		"`-d`",
		"`trace*.txt`",
		"`trace*.csv`",
		"`trace*.tsv`",
		"`.nextflow.log`",
		"`.nextflow_*.log`",
		"`.gosh/index.sqlite`",
		"`gosh index --refresh`",
		"`-with-trace`",
		"`ab/c123def`",
		"`--sample` is a name/tag substring filter",
	}
	for _, want := range required {
		if !strings.Contains(doc, want) {
			t.Fatalf("go-v1-triage.md missing required scope text %q", want)
		}
	}
}

func TestGoV1TriageDocsIncludeCommandExamplesAndMigrationBoundary(t *testing.T) {
	doc := readGoV1TriageDoc(t)

	requiredExamples := []string{
		"gosh status --run-dir /path/to/run",
		"gosh tasks --run-dir /path/to/run --status FAILED",
		"gosh inspect ab/c123def --run-dir /path/to/run",
		"gosh status -d /path/to/run --json",
		"gosh tasks -d /path/to/run --json --status FAILED",
		"gosh inspect ab/c123def -d /path/to/run --json",
		"gosh index --run-dir /path/to/run --refresh",
		"`gosh run`",
		"`gosh debug`",
		"`gosh help`",
		"not v1 commands",
	}
	for _, want := range requiredExamples {
		if !strings.Contains(doc, want) {
			t.Fatalf("go-v1-triage.md missing command example or migration boundary %q", want)
		}
	}
}

func TestGoV1TriageDocsDoNotOverclaimRealRunValidation(t *testing.T) {
	doc := readGoV1TriageDoc(t)

	requiredLimitations := []string{
		"synthetic-only",
		"not real nf-gOS validation",
		"No sanitized real nf-gOS run fixture was provided",
		"timings are local smoke measurements, not benchmark claims",
	}
	for _, want := range requiredLimitations {
		if !strings.Contains(doc, want) {
			t.Fatalf("go-v1-triage.md missing validation limitation %q", want)
		}
	}
}
