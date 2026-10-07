package parser

import (
	"os"
	"path/filepath"
	"testing"
)

func writeTemp(t *testing.T, name, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestParse23andMe(t *testing.T) {
	path := writeTemp(t, "g.txt", "# comment\nrs1\t1\t100\tAA\n\nrs2 2 200 GT\nrs3\t3\t300\t--\n")
	data, err := ParseRawData(path, "23andme")
	if err != nil {
		t.Fatal(err)
	}
	if data["rs1"] != "AA" || data["rs2"] != "GT" || data["rs3"] != "--" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestLivingDNAUses23andMeLayout(t *testing.T) {
	path := writeTemp(t, "g.txt", "rs1\t1\t100\tAC\n")
	data, err := ParseRawData(path, "livingdna")
	if err != nil {
		t.Fatal(err)
	}
	if data["rs1"] != "AC" {
		t.Fatalf("livingdna parse failed: %v", data)
	}
}

func TestParseAncestry(t *testing.T) {
	path := writeTemp(t, "g.txt", "rsid\tchromosome\tposition\tallele1\tallele2\nrs1\t1\t100\tA\tC\n")
	data, err := ParseRawData(path, "ancestry")
	if err != nil {
		t.Fatal(err)
	}
	if data["rs1"] != "AC" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestParseFTDNA(t *testing.T) {
	path := writeTemp(t, "g.csv", "RSID,CHROM,POS,RESULT\n\"rs1\",\"1\",100,\"CT\"\n")
	data, err := ParseRawData(path, "ftdna")
	if err != nil {
		t.Fatal(err)
	}
	if data["rs1"] != "CT" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestParseFTDNA2(t *testing.T) {
	path := writeTemp(t, "g.csv", "rs1,1,100,C,T\n")
	data, err := ParseRawData(path, "ftdna2")
	if err != nil {
		t.Fatal(err)
	}
	if data["rs1"] != "CT" {
		t.Fatalf("unexpected: %v", data)
	}
}

func TestUnsupportedVendor(t *testing.T) {
	path := writeTemp(t, "g.txt", "rs1\t1\t100\tAA\n")
	if _, err := ParseRawData(path, "nope"); err == nil {
		t.Fatal("expected error for unsupported vendor")
	}
}

func TestIsVCFFile(t *testing.T) {
	vcf := writeTemp(t, "x.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\n")
	if !IsVCFFile(vcf) {
		t.Error("plain VCF not detected")
	}
	txt := writeTemp(t, "x.txt", "rs1\t1\t100\tAA\n")
	if IsVCFFile(txt) {
		t.Error("text file misdetected as VCF")
	}
	if missing := filepath.Join(t.TempDir(), "missing.vcf"); IsVCFFile(missing) {
		t.Error("missing file misdetected as VCF")
	}
}

func TestParseFileAutoDetectsVCF(t *testing.T) {
	path := writeTemp(t, "x.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n1\t100\trs1\tA\tG\t10\tPASS\t.\tGT\t0/1\n")
	data, src, err := ParseFile(path, "23andme", nil, BuildAuto)
	if err != nil {
		t.Fatal(err)
	}
	if src.Vendor != "vcf" {
		t.Fatalf("expected vcf source, got %s", src.Vendor)
	}
	if data["rs1"] != "AG" {
		t.Fatalf("unexpected: %v", data)
	}
}
