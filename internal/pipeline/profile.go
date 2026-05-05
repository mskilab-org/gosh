package pipeline

import (
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
)

// Profile is the internal, v1-private pipeline behavior seam.
//
// It is intentionally not an external plugin API: there is no dynamic loading,
// user-installed profile discovery, network access, or public compatibility
// promise. Future in-tree profiles can customize sample matching/extraction,
// process presentation/grouping, and deterministic error enrichment without
// hard-coding pipeline-specific concepts into trace, task, log, or render
// packages.
type Profile interface {
	// MatchSample reports whether a task-like name/tag pair matches a --sample query.
	MatchSample(input SampleMatchInput) bool
	// ProcessView returns behavior-preserving process display and grouping labels.
	ProcessView(input ProcessViewInput) ProcessView
	// EnrichError returns deterministic, profile-specific error enrichment fields.
	EnrichError(input ErrorEnrichmentInput) (ErrorEnrichment, error)
}

// SampleMatchInput is the profile-facing shape for --sample matching. The
// default v1 profile treats Query as a case-insensitive substring of Name or Tag;
// it does not parse biological sample identifiers.
type SampleMatchInput struct {
	Query string
	Name  string
	Tag   string
}

// ProcessViewInput is the profile-facing shape for process presentation.
type ProcessViewInput struct {
	Task domain.Task
}

// ProcessView is the profile-provided process presentation/grouping result.
type ProcessView struct {
	Display string
	Group   string
}

// ErrorEnrichmentInput contains deterministic evidence available to profile
// error-enrichment hooks. Profiles must not perform network access or inspect
// mutable Nextflow run artifacts beyond data already supplied by callers.
type ErrorEnrichmentInput struct {
	Task            *domain.Task
	LogOnlyEvidence *domain.LogOnlyTaskEvidence
	Diagnostics     []domain.Diagnostic
}

// ErrorEnrichment is the additive result of profile error enrichment. Empty
// fields mean no profile-specific enrichment.
type ErrorEnrichment struct {
	Task            *domain.Task
	LogOnlyEvidence *domain.LogOnlyTaskEvidence
	Diagnostics     []domain.Diagnostic
}

// DefaultProfile is the pipeline-agnostic v1 profile. Its implementation must
// preserve current behavior and avoid biological-sample or nf-gOS-specific
// parsing.
type DefaultProfile struct{}

func (DefaultProfile) MatchSample(input SampleMatchInput) bool {
	query := strings.TrimSpace(input.Query)
	if query == "" {
		return true
	}

	needle := strings.ToLower(query)
	return strings.Contains(strings.ToLower(input.Name), needle) ||
		strings.Contains(strings.ToLower(input.Tag), needle)
}

func (DefaultProfile) ProcessView(input ProcessViewInput) ProcessView {
	return ProcessView{Display: input.Task.Process, Group: input.Task.Process}
}

func (DefaultProfile) EnrichError(input ErrorEnrichmentInput) (ErrorEnrichment, error) {
	return ErrorEnrichment{
		Task:            input.Task,
		LogOnlyEvidence: input.LogOnlyEvidence,
		Diagnostics:     input.Diagnostics,
	}, nil
}

func NewDefaultProfile() Profile {
	return DefaultProfile{}
}

func defaultProfile(profile Profile) Profile {
	if profile == nil {
		return NewDefaultProfile()
	}
	return profile
}

func MatchSample(profile Profile, input SampleMatchInput) bool {
	return defaultProfile(profile).MatchSample(input)
}

func ProcessViewForTask(profile Profile, task domain.Task) ProcessView {
	return defaultProfile(profile).ProcessView(ProcessViewInput{Task: task})
}

func EnrichError(profile Profile, input ErrorEnrichmentInput) (ErrorEnrichment, error) {
	return defaultProfile(profile).EnrichError(input)
}
