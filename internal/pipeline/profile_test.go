package pipeline

import (
	"errors"
	"reflect"
	"testing"

	"github.com/mskilab-org/gosh/internal/domain"
)

func TestNewDefaultProfileReturnsDefaultProfile(t *testing.T) {
	got := NewDefaultProfile()
	if got == nil {
		t.Fatal("NewDefaultProfile() returned nil")
	}
	if _, ok := got.(DefaultProfile); !ok {
		t.Fatalf("NewDefaultProfile() type = %T, want pipeline.DefaultProfile", got)
	}
}

type recordingProfile struct {
	sampleResult bool
	sampleCalls  int
	sampleInput  SampleMatchInput

	processResult ProcessView
	processCalls  int
	processInput  ProcessViewInput

	errorResult ErrorEnrichment
	errorErr    error
	errorCalls  int
	errorInput  ErrorEnrichmentInput
}

func (p *recordingProfile) MatchSample(input SampleMatchInput) bool {
	p.sampleCalls++
	p.sampleInput = input
	return p.sampleResult
}

func (p *recordingProfile) ProcessView(input ProcessViewInput) ProcessView {
	p.processCalls++
	p.processInput = input
	return p.processResult
}

func (p *recordingProfile) EnrichError(input ErrorEnrichmentInput) (ErrorEnrichment, error) {
	p.errorCalls++
	p.errorInput = input
	return p.errorResult, p.errorErr
}

func TestMatchSampleDelegatesToConfiguredProfile(t *testing.T) {
	input := SampleMatchInput{Query: "tumor", Name: "ALIGN (tumor-a)", Tag: "tumor-a"}
	profile := &recordingProfile{sampleResult: false}

	if got := MatchSample(profile, input); got {
		t.Fatal("MatchSample returned true, want configured profile decision false")
	}
	if profile.sampleCalls != 1 {
		t.Fatalf("configured profile MatchSample calls = %d, want 1", profile.sampleCalls)
	}
	if profile.sampleInput != input {
		t.Fatalf("configured profile input = %#v, want %#v", profile.sampleInput, input)
	}
}

func TestMatchSampleFallsBackToDefaultProfileWhenNil(t *testing.T) {
	input := SampleMatchInput{Query: "TUMOR", Name: "ALIGN (tumor-a)", Tag: "normal-a"}

	if got := MatchSample(nil, input); !got {
		t.Fatal("MatchSample(nil, input) returned false, want default profile case-insensitive name/tag match")
	}
}

func TestProcessViewForTaskDelegatesToConfiguredProfile(t *testing.T) {
	task := domain.Task{ID: "aa/bb1234", Process: "PIPE:ALIGN", Name: "PIPE:ALIGN (tumor-a)", Tag: "tumor-a"}
	want := ProcessView{Display: "Alignment", Group: "Alignment steps"}
	profile := &recordingProfile{processResult: want}

	got := ProcessViewForTask(profile, task)
	if got != want {
		t.Fatalf("ProcessViewForTask returned %#v, want configured profile result %#v", got, want)
	}
	if profile.processCalls != 1 {
		t.Fatalf("configured profile ProcessView calls = %d, want 1", profile.processCalls)
	}
	if profile.processInput != (ProcessViewInput{Task: task}) {
		t.Fatalf("configured profile input = %#v, want task %#v", profile.processInput, task)
	}
}

func TestProcessViewForTaskFallsBackToDefaultProfileWhenNil(t *testing.T) {
	task := domain.Task{ID: "aa/bb1234", Process: "PIPE:ALIGN", Name: "PIPE:ALIGN (tumor-a)", Tag: "tumor-a"}
	want := ProcessView{Display: "PIPE:ALIGN", Group: "PIPE:ALIGN"}

	got := ProcessViewForTask(nil, task)
	if got != want {
		t.Fatalf("ProcessViewForTask(nil, task) = %#v, want default process view %#v", got, want)
	}
}

func TestProcessViewForTaskFallsBackToDefaultProfileForEmptyProcess(t *testing.T) {
	got := ProcessViewForTask(nil, domain.Task{})
	want := ProcessView{}
	if got != want {
		t.Fatalf("ProcessViewForTask(nil, empty task) = %#v, want %#v", got, want)
	}
}

func TestEnrichErrorDelegatesToConfiguredProfile(t *testing.T) {
	inputTask := &domain.Task{ID: "aa/bb1234", Process: "ALIGN"}
	inputDiagnostic := domain.Diagnostic{Severity: domain.DiagnosticWarning, Code: "input", Message: "input diagnostic"}
	outputTask := &domain.Task{ID: "cc/dd5678", Process: "REPORT"}
	outputDiagnostic := domain.Diagnostic{Severity: domain.DiagnosticInfo, Code: "profile", Message: "profile diagnostic"}
	want := ErrorEnrichment{Task: outputTask, Diagnostics: []domain.Diagnostic{outputDiagnostic}}
	profile := &recordingProfile{errorResult: want}
	input := ErrorEnrichmentInput{Task: inputTask, Diagnostics: []domain.Diagnostic{inputDiagnostic}}

	got, err := EnrichError(profile, input)
	if err != nil {
		t.Fatalf("EnrichError returned error: %v", err)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("EnrichError returned %#v, want configured profile result %#v", got, want)
	}
	if profile.errorCalls != 1 {
		t.Fatalf("configured profile EnrichError calls = %d, want 1", profile.errorCalls)
	}
	if !reflect.DeepEqual(profile.errorInput, input) {
		t.Fatalf("configured profile input = %#v, want %#v", profile.errorInput, input)
	}
}

func TestEnrichErrorFallsBackToDefaultProfileWhenNil(t *testing.T) {
	task := &domain.Task{ID: "aa/bb1234", Process: "ALIGN"}
	logOnlyEvidence := &domain.LogOnlyTaskEvidence{ID: "aa/bb1234", Process: "ALIGN"}
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticWarning, Code: "existing", Message: "kept existing evidence"},
	}

	got, err := EnrichError(nil, ErrorEnrichmentInput{
		Task:            task,
		LogOnlyEvidence: logOnlyEvidence,
		Diagnostics:     diagnostics,
	})
	if err != nil {
		t.Fatalf("EnrichError(nil, input) returned error: %v", err)
	}
	if got.Task != task {
		t.Fatalf("Task pointer = %p, want original %p", got.Task, task)
	}
	if got.LogOnlyEvidence != logOnlyEvidence {
		t.Fatalf("LogOnlyEvidence pointer = %p, want original %p", got.LogOnlyEvidence, logOnlyEvidence)
	}
	if !reflect.DeepEqual(got.Diagnostics, diagnostics) {
		t.Fatalf("Diagnostics = %#v, want %#v", got.Diagnostics, diagnostics)
	}
}

func TestEnrichErrorPropagatesConfiguredProfileError(t *testing.T) {
	wantErr := errors.New("profile enrichment failed")
	profile := &recordingProfile{errorErr: wantErr}

	_, err := EnrichError(profile, ErrorEnrichmentInput{})
	if !errors.Is(err, wantErr) {
		t.Fatalf("EnrichError error = %v, want %v", err, wantErr)
	}
	if profile.errorCalls != 1 {
		t.Fatalf("configured profile EnrichError calls = %d, want 1", profile.errorCalls)
	}
}

func TestNewDefaultProfileProvidesUsableIdentityErrorEnrichment(t *testing.T) {
	task := &domain.Task{ID: "aa/bb1234", Process: "ALIGN"}
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticWarning, Code: "example", Message: "kept existing evidence"},
	}

	got, err := NewDefaultProfile().EnrichError(ErrorEnrichmentInput{
		Task:        task,
		Diagnostics: diagnostics,
	})
	if err != nil {
		t.Fatalf("NewDefaultProfile().EnrichError returned error: %v", err)
	}
	if got.Task != task {
		t.Fatalf("Task pointer = %p, want original %p", got.Task, task)
	}
	if !reflect.DeepEqual(got.Diagnostics, diagnostics) {
		t.Fatalf("Diagnostics = %#v, want %#v", got.Diagnostics, diagnostics)
	}
}

func TestDefaultProfileEnrichErrorPreservesSuppliedEvidenceAndDiagnostics(t *testing.T) {
	exitCode := 137
	task := &domain.Task{
		RowOrder:     7,
		ID:           "aa/bb1234",
		Status:       domain.TaskStatusFailed,
		Process:      "ALIGN",
		Name:         "ALIGN (tumor-a)",
		Tag:          "tumor-a",
		Workdir:      "/runs/work/aa/bb1234",
		Exit:         &exitCode,
		ErrorSummary: "command failed",
	}
	logOnlyEvidence := &domain.LogOnlyTaskEvidence{
		ID:             "aa/bb1234",
		Workdir:        "/runs/work/aa/bb1234",
		Process:        "ALIGN",
		Name:           "ALIGN (tumor-a)",
		ObservedStatus: domain.TaskStatusFailed,
		Exit:           &exitCode,
		ErrorSummary:   "command failed",
		ErrorBlock:     "failed block",
		Sources: []domain.LogOnlyEvidenceSource{
			{Kind: domain.LogOnlyEvidenceSourceLog, Path: "/runs/.nextflow.log", Detail: "selected log"},
		},
		Completeness: domain.LogOnlyEvidencePartial,
	}
	diagnostics := []domain.Diagnostic{
		{Severity: domain.DiagnosticWarning, Code: "paired-log-ambiguous", Message: "kept existing evidence", Detail: "diagnostic detail"},
	}

	got, err := (DefaultProfile{}).EnrichError(ErrorEnrichmentInput{
		Task:            task,
		LogOnlyEvidence: logOnlyEvidence,
		Diagnostics:     diagnostics,
	})
	if err != nil {
		t.Fatalf("DefaultProfile.EnrichError returned error: %v", err)
	}
	if got.Task != task {
		t.Fatalf("Task pointer = %p, want original %p", got.Task, task)
	}
	if got.LogOnlyEvidence != logOnlyEvidence {
		t.Fatalf("LogOnlyEvidence pointer = %p, want original %p", got.LogOnlyEvidence, logOnlyEvidence)
	}
	if !reflect.DeepEqual(got.Diagnostics, diagnostics) {
		t.Fatalf("Diagnostics = %#v, want %#v", got.Diagnostics, diagnostics)
	}
}

func TestDefaultProfileEnrichErrorAcceptsEmptyInput(t *testing.T) {
	got, err := (DefaultProfile{}).EnrichError(ErrorEnrichmentInput{})
	if err != nil {
		t.Fatalf("DefaultProfile.EnrichError(empty) returned error: %v", err)
	}
	if got.Task != nil {
		t.Fatalf("Task = %#v, want nil", got.Task)
	}
	if got.LogOnlyEvidence != nil {
		t.Fatalf("LogOnlyEvidence = %#v, want nil", got.LogOnlyEvidence)
	}
	if len(got.Diagnostics) != 0 {
		t.Fatalf("Diagnostics length = %d, want 0 (%#v)", len(got.Diagnostics), got.Diagnostics)
	}
}
