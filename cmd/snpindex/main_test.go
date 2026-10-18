package main

import "testing"

func TestAsmHitPrefersPrimary(t *testing.T) {
	var h asmHit
	h.set("1", 100, 'A')
	h.set("1_KI270760v1_alt", 50, 'T')
	if h.chrom != "1" || h.pos != 100 || h.ref != 'A' {
		t.Fatalf("alt contig overwrote primary: %+v", h)
	}
	// A later primary mapping still replaces an earlier one.
	h.set("1", 200, 'G')
	if h.chrom != "1" || h.pos != 200 || h.ref != 'G' {
		t.Fatalf("primary did not replace primary: %+v", h)
	}
}

func TestAsmHitAltOnlyIsUnusable(t *testing.T) {
	var row snpRow
	row.hg38.set("4_GL383527v1_alt", 77499, 'T')
	row.hg19.set("4", 156834028, 'T')
	chrom, pos19, pos38, ref, ok := row.fields()
	if !ok || chrom != "4" || pos19 != 156834028 || pos38 != 0 || ref != 'T' {
		t.Fatalf("fields = %s %d %d %c %v", chrom, pos19, pos38, ref, ok)
	}
}

func TestIsPrimaryChrom(t *testing.T) {
	for _, c := range []string{"1", "22", "X", "Y", "M"} {
		if !isPrimaryChrom(c) {
			t.Errorf("%s should be primary", c)
		}
	}
	for _, c := range []string{"", "0", "23", "1_KI270760v1_alt", "MT"} {
		if isPrimaryChrom(c) {
			t.Errorf("%s should not be primary", c)
		}
	}
}
