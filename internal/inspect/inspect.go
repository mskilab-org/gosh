package inspect

import (
	"bufio"
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/mskilab-org/gosh/internal/domain"
)

type SnippetOptions struct {
	MaxBytes int64
	MaxLines int
}

func InventoryCommandFiles(ctx context.Context, workdir string, options SnippetOptions) (domain.CommandFileInventory, error) {
	inventory := domain.CommandFileInventory{Workdir: workdir}
	if ctx == nil {
		return inventory, fmt.Errorf("nil context")
	}
	if err := ctx.Err(); err != nil {
		return inventory, err
	}

	kinds := []domain.CommandFileKind{
		domain.CommandFileShell,
		domain.CommandFileLog,
		domain.CommandFileErr,
		domain.CommandFileOut,
		domain.CommandFileRun,
	}
	inventory.Files = make([]domain.CommandFile, 0, len(kinds))

	for _, kind := range kinds {
		if err := ctx.Err(); err != nil {
			return inventory, err
		}

		path := CommandFilePath(workdir, kind)
		file := domain.CommandFile{Kind: kind, Path: path}
		info, err := os.Stat(path)
		if err != nil {
			if os.IsNotExist(err) {
				inventory.Files = append(inventory.Files, file)
				continue
			}
			return inventory, fmt.Errorf("stat command file %q: %w", path, err)
		}
		if info.IsDir() {
			return inventory, fmt.Errorf("command file path %q is a directory", path)
		}

		file.Exists = true
		file.Size = info.Size()
		snippet, err := ReadBoundedSnippet(ctx, path, SnippetStrategyForKind(kind), options)
		if err != nil {
			return inventory, fmt.Errorf("read command file snippet %q: %w", path, err)
		}
		file.Snippet = snippet
		inventory.Files = append(inventory.Files, file)
	}

	return inventory, nil
}

func CommandFilePath(workdir string, kind domain.CommandFileKind) string {
	return filepath.Join(workdir, string(kind))
}

func ReadBoundedSnippet(ctx context.Context, path string, strategy domain.SnippetStrategy, options SnippetOptions) (*domain.Snippet, error) {
	if ctx == nil {
		return nil, fmt.Errorf("nil context")
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	if options.MaxBytes <= 0 {
		return nil, fmt.Errorf("max bytes must be positive")
	}
	if options.MaxLines <= 0 {
		return nil, fmt.Errorf("max lines must be positive")
	}
	if strategy != domain.SnippetStrategyHead && strategy != domain.SnippetStrategyTail && strategy != domain.SnippetStrategyError {
		return nil, fmt.Errorf("unsupported snippet strategy %q", strategy)
	}

	info, err := os.Stat(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("missing command file %q: %w", path, err)
		}
		return nil, fmt.Errorf("stat command file %q: %w", path, err)
	}
	if info.IsDir() {
		return nil, fmt.Errorf("command file path %q is a directory", path)
	}

	file, err := os.Open(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("missing command file %q: %w", path, err)
		}
		return nil, fmt.Errorf("open command file %q: %w", path, err)
	}
	defer file.Close()

	lines, err := readSnippetLines(ctx, file, options.MaxBytes)
	if err != nil {
		return nil, fmt.Errorf("read command file %q: %w", path, err)
	}

	selected, errorMatched := selectSnippetLines(lines, strategy, options.MaxLines)
	content, included, bytesTruncated := boundedSnippetContent(selected, options.MaxBytes, strategy, errorMatched)
	truncated := bytesTruncated || !coversWholeFile(included, len(lines))

	snippet := &domain.Snippet{
		Path:      path,
		Strategy:  strategy,
		Content:   content,
		Truncated: truncated,
		MaxBytes:  options.MaxBytes,
	}
	if len(included) > 0 {
		snippet.StartLine = included[0].Number
		snippet.EndLine = included[len(included)-1].Number
	}
	return snippet, nil
}

type snippetLine struct {
	Number int
	Text   string
}

func readSnippetLines(ctx context.Context, file *os.File, maxBytes int64) ([]snippetLine, error) {
	scanner := bufio.NewScanner(file)
	scanner.Buffer(make([]byte, 0, 64*1024), scannerTokenLimit(maxBytes))

	var lines []snippetLine
	lineNumber := 0
	for scanner.Scan() {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		lineNumber++
		lines = append(lines, snippetLine{Number: lineNumber, Text: strings.TrimSuffix(scanner.Text(), "\r")})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	return lines, nil
}

func scannerTokenLimit(maxBytes int64) int {
	const defaultLimit = 1024 * 1024
	limit := int64(defaultLimit)
	if maxBytes > limit {
		limit = maxBytes
	}
	maxInt := int64(^uint(0) >> 1)
	if limit > maxInt {
		return int(maxInt)
	}
	return int(limit)
}

func selectSnippetLines(lines []snippetLine, strategy domain.SnippetStrategy, maxLines int) ([]snippetLine, bool) {
	if len(lines) == 0 || maxLines <= 0 {
		return nil, false
	}
	if maxLines > len(lines) {
		maxLines = len(lines)
	}

	switch strategy {
	case domain.SnippetStrategyHead:
		return append([]snippetLine(nil), lines[:maxLines]...), false
	case domain.SnippetStrategyTail:
		return append([]snippetLine(nil), lines[len(lines)-maxLines:]...), false
	case domain.SnippetStrategyError:
		for index, line := range lines {
			if lineLooksErrorFocused(line.Text) {
				before := (maxLines - 1) / 2
				start := index - before
				if start < 0 {
					start = 0
				}
				end := start + maxLines
				if end > len(lines) {
					end = len(lines)
					start = end - maxLines
					if start < 0 {
						start = 0
					}
				}
				return append([]snippetLine(nil), lines[start:end]...), true
			}
		}
		return append([]snippetLine(nil), lines[len(lines)-maxLines:]...), false
	default:
		return nil, false
	}
}

func lineLooksErrorFocused(line string) bool {
	lower := strings.ToLower(line)
	for _, term := range []string{"error", "failed", "failure", "exception", "traceback", "terminated", "exit status", "exit code", "killed", "fatal"} {
		if strings.Contains(lower, term) {
			return true
		}
	}
	return false
}

func boundedSnippetContent(lines []snippetLine, maxBytes int64, strategy domain.SnippetStrategy, errorMatched bool) (string, []snippetLine, bool) {
	if strategy == domain.SnippetStrategyTail || (strategy == domain.SnippetStrategyError && !errorMatched) {
		return boundedSnippetContentFromTail(lines, maxBytes)
	}
	return boundedSnippetContentFromHead(lines, maxBytes)
}

func boundedSnippetContentFromHead(lines []snippetLine, maxBytes int64) (string, []snippetLine, bool) {
	var builder strings.Builder
	included := make([]snippetLine, 0, len(lines))
	usedBytes := int64(0)
	for _, line := range lines {
		separatorBytes := int64(0)
		if len(included) > 0 {
			separatorBytes = 1
		}
		lineBytes := int64(len([]byte(line.Text)))
		if usedBytes+separatorBytes+lineBytes > maxBytes {
			if len(included) == 0 && maxBytes > 0 {
				builder.WriteString(truncateStringBytes(line.Text, maxBytes))
				included = append(included, line)
			}
			return builder.String(), included, true
		}
		if separatorBytes > 0 {
			builder.WriteByte('\n')
		}
		builder.WriteString(line.Text)
		included = append(included, line)
		usedBytes += separatorBytes + lineBytes
	}
	return builder.String(), included, false
}

func boundedSnippetContentFromTail(lines []snippetLine, maxBytes int64) (string, []snippetLine, bool) {
	includedReversed := make([]snippetLine, 0, len(lines))
	usedBytes := int64(0)
	bytesTruncated := false
	for index := len(lines) - 1; index >= 0; index-- {
		line := lines[index]
		separatorBytes := int64(0)
		if len(includedReversed) > 0 {
			separatorBytes = 1
		}
		lineBytes := int64(len([]byte(line.Text)))
		if usedBytes+separatorBytes+lineBytes > maxBytes {
			if len(includedReversed) == 0 && maxBytes > 0 {
				includedReversed = append(includedReversed, snippetLine{Number: line.Number, Text: truncateStringBytesFromEnd(line.Text, maxBytes)})
			}
			bytesTruncated = true
			break
		}
		includedReversed = append(includedReversed, line)
		usedBytes += separatorBytes + lineBytes
	}

	included := make([]snippetLine, len(includedReversed))
	for index := range includedReversed {
		included[len(includedReversed)-1-index] = includedReversed[index]
	}
	parts := make([]string, len(included))
	for index, line := range included {
		parts[index] = line.Text
	}
	return strings.Join(parts, "\n"), included, bytesTruncated
}

func coversWholeFile(included []snippetLine, totalLines int) bool {
	if totalLines == 0 {
		return len(included) == 0
	}
	if len(included) != totalLines || len(included) == 0 {
		return false
	}
	return included[0].Number == 1 && included[len(included)-1].Number == totalLines
}

func truncateStringBytes(value string, maxBytes int64) string {
	if int64(len([]byte(value))) <= maxBytes {
		return value
	}
	return string([]byte(value)[:maxBytes])
}

func truncateStringBytesFromEnd(value string, maxBytes int64) string {
	if int64(len([]byte(value))) <= maxBytes {
		return value
	}
	bytes := []byte(value)
	return string(bytes[int64(len(bytes))-maxBytes:])
}

func SnippetStrategyForKind(kind domain.CommandFileKind) domain.SnippetStrategy {
	switch kind {
	case domain.CommandFileShell, domain.CommandFileRun:
		return domain.SnippetStrategyHead
	case domain.CommandFileLog, domain.CommandFileErr:
		return domain.SnippetStrategyError
	case domain.CommandFileOut:
		return domain.SnippetStrategyTail
	default:
		return domain.SnippetStrategyTail
	}
}
