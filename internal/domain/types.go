package domain

import (
	"errors"
	"time"
)

var ErrNotImplemented = errors.New("not implemented")

type RunDir struct {
	Path string
}

type ResultsDir struct {
	Path string
}

type IndexMode string

const (
	IndexModeTraceBacked IndexMode = "trace-backed"
	IndexModeLogOnly     IndexMode = "log-only"
	IndexModeUnsupported IndexMode = "unsupported"
)

type SourceKind string

const (
	SourceKindTrace SourceKind = "trace"
	SourceKindLog   SourceKind = "log"
)

type ArtifactSearchLocation struct {
	Kind        SourceKind
	BaseDir     string
	Patterns    []string
	Description string
}

type SourceFingerprint struct {
	Kind    SourceKind
	Path    string
	ModTime time.Time
	Size    int64
}

type ArtifactSet struct {
	RunDir           RunDir
	ResultsDir       ResultsDir
	Mode             IndexMode
	Trace            *SourceFingerprint
	Log              *SourceFingerprint
	SelectedAt       time.Time
	SearchedPatterns []string
	SearchLocations  []ArtifactSearchLocation
	Diagnostics      []Diagnostic
}

type DiagnosticSeverity string

const (
	DiagnosticInfo    DiagnosticSeverity = "info"
	DiagnosticWarning DiagnosticSeverity = "warning"
	DiagnosticError   DiagnosticSeverity = "error"
)

type Diagnostic struct {
	Severity DiagnosticSeverity
	Code     string
	Message  string
	Detail   string
}

type DiagnosticContextLine struct {
	Label string
	Value string
}

type DiagnosticBlock struct {
	Severity DiagnosticSeverity
	Code     string
	Title    string
	Context  []DiagnosticContextLine
	Details  []string
	Hints    []string
}

type IndexFreshness string

const (
	IndexFreshnessUnknown     IndexFreshness = "unknown"
	IndexFreshnessMissing     IndexFreshness = "missing"
	IndexFreshnessFresh       IndexFreshness = "fresh"
	IndexFreshnessStale       IndexFreshness = "stale"
	IndexFreshnessUnsupported IndexFreshness = "unsupported"
)

type IndexMetadata struct {
	SchemaVersion int
	RunDir        string
	IndexPath     string
	Mode          IndexMode
	Trace         *SourceFingerprint
	Log           *SourceFingerprint
	BuiltAt       time.Time
	Freshness     IndexFreshness
	StaleReason   string
	TaskCount     int
}

type TaskStatus string

const (
	TaskStatusUnknown   TaskStatus = "UNKNOWN"
	TaskStatusSubmitted TaskStatus = "SUBMITTED"
	TaskStatusRunning   TaskStatus = "RUNNING"
	TaskStatusCompleted TaskStatus = "COMPLETED"
	TaskStatusCached    TaskStatus = "CACHED"
	TaskStatusFailed    TaskStatus = "FAILED"
	TaskStatusAborted   TaskStatus = "ABORTED"
)

type Task struct {
	RowOrder int64
	ID       string
	Status   TaskStatus
	Process  string
	Name     string
	Tag      string
	// Workdir is empty when unknown, otherwise a path to the task work directory.
	// When resolved from a short Nextflow hash prefix, it may point to a longer
	// existing on-disk directory while ID remains the canonical trace/log hash.
	Workdir      string
	Exit         *int
	Duration     string
	Realtime     string
	CPUs         string
	Memory       string
	ErrorSummary string
}

type TaskQuery struct {
	ProcessSubstring string
	NameSubstring    string
	SampleSubstring  string
	Status           TaskStatus
	StatusRaw        string
}

type StatusCount struct {
	Status TaskStatus
	Count  int
}

type FailedTaskPreview struct {
	ID           string
	Status       TaskStatus
	Process      string
	Name         string
	Tag          string
	Workdir      string
	Exit         *int
	ErrorSummary string
}

type StatusSummary struct {
	RunDir          RunDir
	Mode            IndexMode
	IndexPath       string
	Freshness       IndexFreshness
	BuiltAt         *time.Time
	Sources         ArtifactSet
	Counts          []StatusCount
	FailedCount     int
	FailedPreview   []FailedTaskPreview
	LogOnlyFailures []LogOnlyFailure
	LogOnlyEvidence []LogOnlyTaskEvidence
	Diagnostics     []Diagnostic
}

type SelectorResolutionKind string

const (
	SelectorResolutionExact     SelectorResolutionKind = "exact"
	SelectorResolutionAmbiguous SelectorResolutionKind = "ambiguous"
	SelectorResolutionNotFound  SelectorResolutionKind = "not-found"
)

type SelectorResolution struct {
	Kind        SelectorResolutionKind
	Selector    string
	Task        *Task
	Matches     []Task
	Diagnostics []Diagnostic
}

type CommandFileKind string

const (
	CommandFileShell CommandFileKind = ".command.sh"
	CommandFileLog   CommandFileKind = ".command.log"
	CommandFileErr   CommandFileKind = ".command.err"
	CommandFileOut   CommandFileKind = ".command.out"
	CommandFileRun   CommandFileKind = ".command.run"
)

type SnippetStrategy string

const (
	SnippetStrategyHead  SnippetStrategy = "head"
	SnippetStrategyTail  SnippetStrategy = "tail"
	SnippetStrategyError SnippetStrategy = "error-focused"
)

type Snippet struct {
	Path      string
	Strategy  SnippetStrategy
	StartLine int
	EndLine   int
	Content   string
	Truncated bool
	MaxBytes  int64
}

type CommandFile struct {
	Kind    CommandFileKind
	Path    string
	Exists  bool
	Size    int64
	Snippet *Snippet
}

type CommandFileInventory struct {
	Workdir string
	Files   []CommandFile
}

type TaskDossier struct {
	Task        Task
	Inventory   CommandFileInventory
	Diagnostics []Diagnostic
}

type LogOnlyFailure struct {
	ID string
	// Workdir is empty when unknown, otherwise a path to the task work directory.
	// When resolved from a short Nextflow hash prefix, it may point to a longer
	// existing on-disk directory while ID remains the canonical trace/log hash.
	Workdir      string
	Process      string
	Name         string
	Exit         *int
	ErrorSummary string
	ErrorBlock   string
}

type LogOnlyEvidenceSourceKind string

const (
	LogOnlyEvidenceSourceLog     LogOnlyEvidenceSourceKind = "log"
	LogOnlyEvidenceSourceWorkdir LogOnlyEvidenceSourceKind = "workdir"
	LogOnlyEvidenceSourceCommand LogOnlyEvidenceSourceKind = "command-file"
)

type LogOnlyEvidenceCompleteness string

const (
	LogOnlyEvidencePartial LogOnlyEvidenceCompleteness = "partial"
	LogOnlyEvidenceNone    LogOnlyEvidenceCompleteness = "none"
)

type LogOnlyEvidenceSource struct {
	Kind   LogOnlyEvidenceSourceKind
	Path   string
	Detail string
}

type LogOnlyTaskEvidence struct {
	ID string
	// Workdir is empty when unknown, otherwise a path to the task work directory.
	// When resolved from a short Nextflow hash prefix, it may point to a longer
	// existing on-disk directory while ID remains the canonical trace/log hash.
	Workdir               string
	Process               string
	Name                  string
	ObservedStatus        TaskStatus
	Exit                  *int
	ErrorSummary          string
	ErrorBlock            string
	Sources               []LogOnlyEvidenceSource
	Completeness          LogOnlyEvidenceCompleteness
	CommandFilesAvailable bool
}

type IndexDiagnostics struct {
	RunDir      RunDir
	Artifacts   ArtifactSet
	Metadata    *IndexMetadata
	Diagnostics []Diagnostic
}

type OutputFormat string

const (
	OutputFormatHuman OutputFormat = "human"
	OutputFormatJSON  OutputFormat = "json"
)

type StatusView struct {
	Summary StatusSummary
	Format  OutputFormat
}

type TasksView struct {
	Tasks       []Task
	Query       TaskQuery
	Metadata    *IndexMetadata
	Diagnostics []Diagnostic
	Format      OutputFormat
}

type InspectEvidenceKind string

const (
	InspectEvidenceTraceBacked InspectEvidenceKind = "trace-backed"
	InspectEvidenceLogOnly     InspectEvidenceKind = "log-only-partial"
)

type LogOnlySelectorResolution struct {
	Kind        SelectorResolutionKind
	Selector    string
	Evidence    *LogOnlyTaskEvidence
	Matches     []LogOnlyTaskEvidence
	Diagnostics []Diagnostic
}

type LogOnlyTaskDossier struct {
	Evidence    LogOnlyTaskEvidence
	Inventory   CommandFileInventory
	Diagnostics []Diagnostic
}

type InspectView struct {
	Resolution        SelectorResolution
	Dossier           *TaskDossier
	EvidenceKind      InspectEvidenceKind
	LogOnlyResolution *LogOnlySelectorResolution
	LogOnlyDossier    *LogOnlyTaskDossier
	Diagnostics       []Diagnostic
	Format            OutputFormat
}

type IndexView struct {
	Diagnostics IndexDiagnostics
	Format      OutputFormat
}
