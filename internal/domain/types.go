package domain

import (
	"errors"
	"time"
)

var ErrNotImplemented = errors.New("not implemented")

type RunDir struct {
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

type SourceFingerprint struct {
	Kind    SourceKind
	Path    string
	ModTime time.Time
	Size    int64
}

type ArtifactSet struct {
	RunDir           RunDir
	Mode             IndexMode
	Trace            *SourceFingerprint
	Log              *SourceFingerprint
	SelectedAt       time.Time
	SearchedPatterns []string
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
	RowOrder     int64
	ID           string
	Status       TaskStatus
	Process      string
	Name         string
	Tag          string
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
	ID           string
	Workdir      string
	Process      string
	Name         string
	Exit         *int
	ErrorSummary string
	ErrorBlock   string
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

type InspectView struct {
	Resolution  SelectorResolution
	Dossier     *TaskDossier
	Diagnostics []Diagnostic
	Format      OutputFormat
}

type IndexView struct {
	Diagnostics IndexDiagnostics
	Format      OutputFormat
}
