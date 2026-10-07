package parser

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"

	"github.com/emre-tarhan/admixture/pkg/models"
)

// wantedFor builds a Wanted set matching the synthetic fixtures below.
func wantedFor(t testing.TB) *models.Wanted {
	t.Helper()
	m := &models.Model{
		SNPs: []string{"rs1", "rs2", "rsPosOnly"},
		Coords: []*models.Coord{
			{RsID: "rs1", Chrom: "1", Pos38: 100, Pos19: 150, Ref: 'A'},
			{RsID: "rs2", Chrom: "1", Pos38: 200, Pos19: 250, Ref: 'C'},
			{RsID: "rsPosOnly", Chrom: "2", Pos38: 5000, Pos19: 5100, Ref: 'G'},
		},
	}
	return models.BuildWanted([]*models.Model{m})
}

const vcfHeader = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=chr1,length=248956422>\n" + // hg38 chr1
	"##contig=<ID=chr2,length=242193529>\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

func TestParseVCFRsIDMatching(t *testing.T) {
	body := vcfHeader +
		"chr1\t100\trs1\tA\tG\t10\tPASS\t.\tGT\t0/1\n" +
		"chr1\t200\trs2\tC\tA,T\t10\tPASS\t.\tGT\t1/2\n" +
		"chr1\t300\trs3\tT\tC\t10\tPASS\t.\tGT\t1|0\n"
	data, build, stats, err := ParseVCF(writeVCF(t, body), nil, BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if build != Build38 {
		t.Fatalf("build = %q, want hg38", build)
	}
	if data["rs1"] != "AG" || data["rs2"] != "AT" || data["rs3"] != "CT" {
		t.Fatalf("unexpected: %v", data)
	}
	if stats.Kept != 3 {
		t.Fatalf("kept = %d, want 3", stats.Kept)
	}
}

func TestParseVCFPositionMatching(t *testing.T) {
	// No rsIDs at all (like WGS callsets): position keys must appear.
	body := vcfHeader +
		"chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT:GQ\t1/1:30\n" +
		"chr2\t5000\t.\tG\tT\t10\tPASS\t.\tGT:GQ\t0/1:30\n"
	data, build, _, err := ParseVCF(writeVCF(t, body), wantedFor(t), BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if build != Build38 {
		t.Fatalf("build = %q, want hg38", build)
	}
	if data["hg38:1:100"] != "GG" {
		t.Fatalf("missing hg38 position key: %v", data)
	}
	if data["hg38:2:5000"] != "GT" {
		t.Fatalf("missing hg38 position key: %v", data)
	}
}

func TestParseVCFBuildForced(t *testing.T) {
	// Header without contigs: build unknown, forced hg19 must be used.
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n1\t150\t.\tA\tG\t10\tPASS\t.\tGT\t0/1\n"
	data, build, _, err := ParseVCF(writeVCF(t, body), wantedFor(t), Build19)
	if err != nil {
		t.Fatal(err)
	}
	if build != Build19 {
		t.Fatalf("build = %q, want hg19", build)
	}
	if data["hg19:1:150"] != "AG" {
		t.Fatalf("missing hg19 position key: %v", data)
	}
}

func TestParseVCFBuildUndetermined(t *testing.T) {
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n1\t100\trs1\tA\tG\t10\tPASS\t.\tGT\t0/1\n"
	data, build, _, err := ParseVCF(writeVCF(t, body), nil, BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if build != "" {
		t.Fatalf("build = %q, want empty", build)
	}
	// rsID matching still works.
	if data["rs1"] != "AG" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestParseVCFGVCFCRefBlocks(t *testing.T) {
	// Ref block spanning 4000..6000 on chr2 covers rsPosOnly (hg38 5000).
	body := vcfHeader +
		"chr2\t4000\t.\tG\t<*>\t0\t.\tEND=6000\tGT:GQ\t0/0:20\n" +
		"chr2\t4100\t.\tT\t<*>\t0\t.\tEND=4200\tGT:GQ\t./.:20\n" // no-call block
	data, _, stats, err := ParseVCF(writeVCF(t, body), wantedFor(t), BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if data["hg38:2:5000"] != "GG" {
		t.Fatalf("ref block did not resolve hom-ref SNP: %v", data)
	}
	if stats.BlockCoords != 1 {
		t.Fatalf("blockCoords = %d, want 1", stats.BlockCoords)
	}
}

func TestParseVCFSkips(t *testing.T) {
	body := vcfHeader +
		"chr1\t110\t.\tA\tG\t10\tLowQual\t.\tGT\t1/1\n" + // filtered
		"chr1\t120\t.\tAT\tG\t10\tPASS\t.\tGT\t1/1\n" + // indel
		"chr1\t130\t.\tA\t<DEL>\t10\tPASS\t.\tGT\t1/1\n" + // symbolic
		"chr1\t140\t.\tA\tG\t10\tPASS\t.\tGT\t0/.\n" + // partial no-call
		"chr1\t150\t.\tA\tG\t10\tPASS\t.\tGL\t0.1,0.2\n" + // no GT field
		"chr1\t160\t.\tA\tG,*\t10\tPASS\t.\tGT\t1/2\n" // spanning deletion allele
	data, _, _, err := ParseVCF(writeVCF(t, body), nil, BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if len(data) != 0 {
		t.Fatalf("expected nothing parsed, got %v", data)
	}
}

func TestParseVCFHemizygous(t *testing.T) {
	// Haploid GT (chrX male): single allele is kept as-is.
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchrX\t900\trsX\tC\tA\t10\tPASS\t.\tGT\t1\n"
	data, _, _, err := ParseVCF(writeVCF(t, body), nil, BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if data["rsX"] != "A" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestParseVCFGzipped(t *testing.T) {
	body := vcfHeader + "chr1\t100\trs1\tA\tG\t10\tPASS\t.\tGT\t1/1\n"
	path := filepath.Join(t.TempDir(), "s.vcf.gz")
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write([]byte(body)); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, buf.Bytes(), 0644); err != nil {
		t.Fatal(err)
	}

	data, _, _, err := ParseVCF(path, nil, BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if data["rs1"] != "GG" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestParseVCFNormalizeChrom(t *testing.T) {
	if got := normalizeChrom("chr23"); got != "X" {
		t.Fatalf("chr23 -> %q", got)
	}
	if got := normalizeChrom("MT"); got != "M" {
		t.Fatalf("MT -> %q", got)
	}
	if got := normalizeChrom("chr12"); got != "12" {
		t.Fatalf("chr12 -> %q", got)
	}
}

func TestResolveGenotype(t *testing.T) {
	cases := []struct {
		gt, ref, alt string
		want         string
		ok           bool
	}{
		{"0/0", "A", "G", "AA", true},
		{"0|1", "A", "G", "AG", true},
		{"1/1", "A", "G", "GG", true},
		{"1/2", "A", "G,T", "GT", true},
		{"./.", "A", "G", "", false},
		{"0/.", "A", "G", "", false},
		{".", "A", "G", "", false},
		{"1/2", "A", "G,<NON_REF>", "", false}, // symbolic allele in GT
	}
	for _, c := range cases {
		got, ok := resolveGenotype(c.gt, c.ref, c.alt)
		if ok != c.ok || (ok && got != c.want) {
			t.Errorf("resolveGenotype(%q,%q,%q) = %q,%v; want %q,%v", c.gt, c.ref, c.alt, got, ok, c.want, c.ok)
		}
	}
}

func writeVCF(t *testing.T, body string) string {
	t.Helper()
	return writeTemp(t, "in.vcf", body)
}
