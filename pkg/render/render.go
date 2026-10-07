// Package render draws human-friendly terminal output for admixture
// results — bar charts instead of bare numbers. It uses only ANSI escape
// codes, no external dependencies, and degrades to plain text when color
// is unavailable.
package render

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/emre-tarhan/admixture/pkg/result"
)

const (
	barWidth    = 24
	minSharePct = 0.01 // hide populations below 0.01%, like the text output
)

// Pretty writes results as bar charts sorted by share, descending.
func Pretty(w io.Writer, results []result.Result) {
	color := SupportsColor(os.Stdout)
	for i, r := range results {
		if i > 0 {
			fmt.Fprintln(w)
		}
		writeModel(w, r, color)
	}
}

func writeModel(w io.Writer, r result.Result, color bool) {
	title := fmt.Sprintf(" %s · %s / %s SNPs matched", r.Model, formatInt(r.MatchedSNPs), formatInt(r.TotalSNPs))
	if color {
		title = "\x1b[1m" + title + "\x1b[0m"
	}
	fmt.Fprintln(w, title)
	fmt.Fprintln(w, strings.Repeat("─", 64))

	entries := r.Sorted(minSharePct / 100)
	if len(entries) == 0 {
		fmt.Fprintln(w, "  (no populations above 0.01%)")
		return
	}

	max := entries[0].Share
	nameWidth := 0
	for _, e := range entries {
		if l := len(e.Name); l > nameWidth {
			nameWidth = l
		}
	}

	for _, e := range entries {
		pct := e.Share * 100
		b := bar(e.Share/max, barWidth, color)
		fmt.Fprintf(w, "  %-*s %s %6.2f%%\n", nameWidth, e.Name, b, pct)
	}
}

// bar draws a relative bar: filled portion proportional to fraction (0..1),
// remainder as light shade.
func bar(fraction float64, width int, color bool) string {
	if fraction < 0 {
		fraction = 0
	}
	if fraction > 1 {
		fraction = 1
	}
	filled := int(fraction*float64(width) + 0.5)
	if filled == 0 && fraction > 0 {
		filled = 1
	}
	b := strings.Repeat("█", filled) + strings.Repeat("░", width-filled)
	if color {
		return "\x1b[32m" + b + "\x1b[0m"
	}
	return b
}

// SupportsColor reports whether ANSI colors should be used on f: it must
// be a character device and the environment must not disable them.
func SupportsColor(f *os.File) bool {
	if os.Getenv("NO_COLOR") != "" {
		return false
	}
	if os.Getenv("TERM") == "dumb" {
		return false
	}
	info, err := f.Stat()
	if err != nil {
		return false
	}
	return info.Mode()&os.ModeCharDevice != 0
}

// formatInt renders n with thousands separators (121334 -> "121,334").
func formatInt(n int) string {
	s := fmt.Sprint(n)
	if len(s) <= 3 {
		return s
	}
	var b strings.Builder
	pre := len(s) % 3
	if pre > 0 {
		b.WriteString(s[:pre])
	}
	for i := pre; i < len(s); i += 3 {
		if b.Len() > 0 {
			b.WriteByte(',')
		}
		b.WriteString(s[i : i+3])
	}
	return b.String()
}
