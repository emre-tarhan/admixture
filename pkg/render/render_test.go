package render

import (
	"bytes"
	"strings"
	"testing"

	"github.com/emre-tarhan/admixture/pkg/result"
)

func TestFormatInt(t *testing.T) {
	cases := map[int]string{
		0:         "0",
		7:         "7",
		123:       "123",
		1234:      "1,234",
		121334:    "121,334",
		1234567:   "1,234,567",
		123456789: "123,456,789",
	}
	for in, want := range cases {
		if got := formatInt(in); got != want {
			t.Errorf("formatInt(%d) = %q, want %q", in, got, want)
		}
	}
}

func TestBar(t *testing.T) {
	if got := bar(0.5, 10, false); got != "█████░░░░░" {
		t.Errorf("bar(0.5) = %q", got)
	}
	if got := bar(1.0, 4, false); got != "████" {
		t.Errorf("bar(1.0) = %q", got)
	}
	if got := bar(0.001, 24, false); !strings.HasPrefix(got, "█") {
		t.Errorf("small positive share should show at least one block, got %q", got)
	}
	if got := bar(2.0, 4, false); got != "████" {
		t.Errorf("bar clamps above 1: %q", got)
	}
}

func TestPrettyOutput(t *testing.T) {
	res := []result.Result{
		{
			Model:       "demo",
			Proportions: map[string]float64{"Pop A": 0.5, "Pop B": 0.3, "Tiny": 0.00001},
			MatchedSNPs: 1234,
			TotalSNPs:   2000,
		},
	}
	var buf bytes.Buffer
	Pretty(&buf, res)
	out := buf.String()

	if !strings.Contains(out, "demo · 1,234 / 2,000 SNPs matched") {
		t.Errorf("header missing: %q", out)
	}
	if !strings.Contains(out, "Pop A") || !strings.Contains(out, "50.00%") {
		t.Errorf("Pop A row missing: %q", out)
	}
	if strings.Contains(out, "Tiny") {
		t.Error("sub-threshold population should be hidden")
	}
	if idxA, idxB := strings.Index(out, "Pop A"), strings.Index(out, "Pop B"); idxA > idxB {
		t.Error("Pop A (larger) should be listed first")
	}
	// Color is disabled in tests (stdout is not a char device in-process).
	if strings.Contains(out, "\x1b[") {
		t.Error("unexpected ANSI colors")
	}
}

func TestPrettyEmpty(t *testing.T) {
	var buf bytes.Buffer
	Pretty(&buf, nil)
	if strings.TrimSpace(buf.String()) != "" {
		t.Errorf("expected empty output, got %q", buf.String())
	}
}
