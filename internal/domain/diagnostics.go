package domain

func NextflowTraceRecommendationDiagnostic() Diagnostic {
	return Diagnostic{
		Severity: DiagnosticInfo,
		Code:     "nextflow_with_trace_recommended",
		Message:  "Run future Nextflow workflows with -with-trace",
		Detail:   "Use `nextflow run ... -with-trace` for future runs so gosh can build a complete trace-backed task index.",
	}
}
