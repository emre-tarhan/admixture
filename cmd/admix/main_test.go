package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// tinyGenome builds a small 23andMe-format file from real K12b rsIDs so the
// integration test exercises the full CLI pipeline.
func tinyGenome(t *testing.T, n int) string {
	t.Helper()
	alleles, err := os.ReadFile("../../pkg/models/data/K12b.alleles")
	if err != nil {
		t.Skip("model data not available")
	}
	var b strings.Builder
	count := 0
	for _, line := range strings.Split(string(alleles), "\n") {
		parts := strings.Fields(line)
		if len(parts) < 3 {
			continue
		}
		b.WriteString(parts[0] + "\t1\t100\t" + parts[2] + parts[2] + "\n") // hom-major
		if count++; count >= n {
			break
		}
	}
	path := filepath.Join(t.TempDir(), "genome.txt")
	if err := os.WriteFile(path, []byte(b.String()), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestRunTextOutput(t *testing.T) {
	path := tinyGenome(t, 200)
	out := capture(t, func() error {
		return run([]string{"-f", path, "-m", "K12b", "-v", "23andme"})
	})
	if !strings.Contains(out, "Model: K12b") {
		t.Fatalf("model header missing: %s", out)
	}
	if !strings.Contains(out, "SNPs matched)") {
		t.Fatalf("stats line missing: %s", out)
	}
}

func TestRunPrettyOutput(t *testing.T) {
	path := tinyGenome(t, 200)
	out := capture(t, func() error {
		return run([]string{"-f", path, "-m", "K12b", "-v", "23andme", "-pretty"})
	})
	if !strings.Contains(out, "K12b ·") {
		t.Fatalf("pretty header missing: %s", out)
	}
	if !strings.Contains(out, "█") {
		t.Fatalf("bar chart missing: %s", out)
	}
}

func TestRunJSONOutput(t *testing.T) {
	path := tinyGenome(t, 100)
	outPath := filepath.Join(t.TempDir(), "out.json")
	if err := run([]string{"-f", path, "-m", "K12b", "--json", "-o", outPath}); err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatal(err)
	}
	s := string(b)
	if !strings.Contains(s, `"model": "K12b"`) || !strings.Contains(s, `"matched_snps"`) || !strings.Contains(s, `"proportions"`) {
		t.Fatalf("JSON shape unexpected: %s", s[:min(400, len(s))])
	}
}

func TestRunMutuallyExclusiveFlags(t *testing.T) {
	path := tinyGenome(t, 10)
	err := run([]string{"-f", path, "-m", "K12b", "--json", "-pretty"})
	if err == nil || !strings.Contains(err.Error(), "mutually exclusive") {
		t.Fatalf("expected mutual exclusion error, got %v", err)
	}
}

func TestRunMissingFileFlag(t *testing.T) {
	if err := run([]string{}); err == nil {
		t.Fatal("expected error without -f")
	}
}

// capture redirects stdout for the duration of fn. The temp file must be
// closed before the test ends — open handles make Windows TempDir cleanup
// fail with "being used by another process".
func capture(t *testing.T, fn func() error) string {
	t.Helper()
	tmp, err := os.CreateTemp(t.TempDir(), "stdout")
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stdout
	os.Stdout = tmp
	runErr := fn()
	os.Stdout = old
	if err := tmp.Close(); err != nil {
		t.Fatal(err)
	}
	if runErr != nil {
		t.Fatalf("run failed: %v", runErr)
	}
	b, err := os.ReadFile(tmp.Name())
	if err != nil {
		t.Fatal(err)
	}
	return string(b)
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
