package parser

import (
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// benchFile writes a synthetic input file and returns its path.
func benchFile(t testing.TB, name, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

func bench23andMe(n int) string {
	var b strings.Builder
	for i := 0; i < n; i++ {
		fmt.Fprintf(&b, "rs%d\t1\t%d\t%s\n", i, i+100, "AG")
	}
	return b.String()
}

func benchVCF(n int) string {
	var b strings.Builder
	b.WriteString(vcfHeader)
	for i := 0; i < n; i++ {
		fmt.Fprintf(&b, "chr1\t%d\trs%d\tA\tG\t50\tPASS\t.\tGT:DP\t%d:30\n", i+100, i, i%3)
	}
	return b.String()
}

func BenchmarkParse23andMe100k(b *testing.B) {
	path := benchFile(b, "bench.txt", bench23andMe(100000))
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := ParseRawData(path, "23andme"); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkParseVCF100k(b *testing.B) {
	path := benchFile(b, "bench.vcf", benchVCF(100000))
	w := wantedFor(&testing.T{})
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, _, err := ParseVCF(path, w, BuildAuto); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkParseVCFGzip100k(b *testing.B) {
	raw := benchVCF(100000)
	path := filepath.Join(b.TempDir(), "bench.vcf.gz")
	f, _ := os.Create(path)
	zw := gzip.NewWriter(f)
	if _, err := zw.Write([]byte(raw)); err != nil {
		b.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		b.Fatal(err)
	}
	if err := f.Close(); err != nil {
		b.Fatal(err)
	}

	w := wantedFor(&testing.T{})
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, _, err := ParseVCF(path, w, BuildAuto); err != nil {
			b.Fatal(err)
		}
	}
}
