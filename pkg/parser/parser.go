package parser

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/emre-tarhan/admixture/pkg/models"
)

// RefMode controls homozygous-reference calls for model SNPs that a
// variant-only VCF simply does not list.
type RefMode string

const (
	// RefAuto assumes homozygous-reference on a variant-only whole-genome
	// VCF and leaves every other file untouched.
	RefAuto RefMode = "auto"
	// RefOn assumes homozygous-reference for any variant-only VCF, including
	// exomes and single-chromosome slices.
	RefOn RefMode = "on"
	// RefOff never invents a genotype for a site the file does not contain.
	RefOff RefMode = "off"
)

// ParseRefMode parses the -homref / homref values.
func ParseRefMode(s string) (RefMode, error) {
	switch RefMode(s) {
	case "", RefAuto:
		return RefAuto, nil
	case RefOn, RefOff:
		return RefMode(s), nil
	default:
		return "", fmt.Errorf("invalid hom-ref mode %q (want auto, on, or off)", s)
	}
}

// Source describes how an input file was interpreted.
type Source struct {
	Vendor string // "23andme", "vcf", ...
	Build  string // reference build used for position matching (VCF only)
	Stats  VCFStats
}

// ParseFile parses a raw genome data file, auto-detecting VCF/gVCF input
// (by content sniffing) and falling back to the named vendor's text
// format. wanted is the key set to extract from VCF input (rsIDs and
// coordinates for the selected models); build forces the reference build
// ("auto" detects it from contig lengths).
func ParseFile(filename, vendor string, wanted *models.Wanted, build Build, ref RefMode) (map[string]string, Source, error) {
	if IsVCFFile(filename) {
		data, b, stats, err := ParseVCF(filename, wanted, build, ref)
		src := Source{Vendor: "vcf", Build: string(b), Stats: stats}
		return data, src, err
	}
	data, err := parseVendorText(filename, vendor)
	return data, Source{Vendor: vendor}, err
}

// IsVCFFile reports whether the file looks like a (possibly gzipped) VCF
// by checking its first line for the fileformat header.
func IsVCFFile(filename string) bool {
	f, err := os.Open(filename)
	if err != nil {
		return false
	}
	defer func() { _ = f.Close() }()

	br := bufio.NewReader(f)
	var r io.Reader = br
	if magic, err := br.Peek(2); err == nil && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(br)
		if err != nil {
			return false
		}
		defer func() { _ = gz.Close() }()
		r = gz
	}
	scanner := bufio.NewScanner(r)
	if !scanner.Scan() {
		return false
	}
	return strings.HasPrefix(scanner.Text(), "##fileformat=VCF")
}

// ParseRawData parses a raw genome data file in one of the vendor text
// formats. For VCF/gVCF input and format auto-detection use ParseFile.
func ParseRawData(filename, vendor string) (map[string]string, error) {
	return parseVendorText(filename, vendor)
}

func parseVendorText(filename, vendor string) (map[string]string, error) {
	switch vendor {
	case "23andme", "wegene", "livingdna":
		// LivingDNA raw downloads use the same 4-column layout as 23andMe.
		return parse23andMe(filename)
	case "ancestry":
		return parseAncestry(filename)
	case "ftdna", "myheritage":
		return parseFTDNA(filename)
	case "ftdna2":
		return parseFTDNA2(filename)
	case "vcf", "gvcf":
		data, _, _, err := ParseVCF(filename, nil, BuildAuto, RefAuto)
		return data, err
	default:
		return nil, fmt.Errorf("unsupported vendor: %s", vendor)
	}
}

// parse23andMe parses the 23andMe layout: rsid chrom pos genotype,
// whitespace-separated, '#' comments allowed.
func parse23andMe(filename string) (map[string]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() { _ = f.Close() }()

	data := make(map[string]string)
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 64<<10), 4<<20)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}
		// Use Fields to handle tabs or spaces
		parts := strings.Fields(line)
		if len(parts) == 4 {
			// rsid chrom pos genotype
			genotype := parts[3]
			data[parts[0]] = genotype
		}
	}
	return data, scanner.Err()
}

// parseAncestry parses the AncestryDNA layout: rsid chrom pos allele1
// allele2, whitespace-separated, one header line.
func parseAncestry(filename string) (map[string]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() { _ = f.Close() }()

	data := make(map[string]string)
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 64<<10), 4<<20)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}
		if strings.HasPrefix(line, "rsid") { // Header
			continue
		}
		parts := strings.Fields(line)
		if len(parts) == 5 {
			// rsid chrom pos allele1 allele2
			genotype := parts[3] + parts[4]
			data[parts[0]] = genotype
		}
	}
	return data, scanner.Err()
}

// parseFTDNA parses the FamilyTreeDNA/MyHeritage CSV layout:
// RSID,CHROM,POS,RESULT — falling back to whitespace separation, with
// quoted fields tolerated.
func parseFTDNA(filename string) (map[string]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() { _ = f.Close() }()

	data := make(map[string]string)
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 64<<10), 4<<20)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") || strings.HasPrefix(line, "RSID") || line == "" {
			continue
		}

		// Try comma first
		parts := strings.Split(line, ",")
		if len(parts) < 4 {
			// Try whitespace
			parts = strings.Fields(line)
		}

		if len(parts) == 4 {
			// RSID, CHROM, POS, RESULT
			// Remove quotes if any
			rsid := strings.Trim(parts[0], "\"")
			geno := strings.Trim(parts[3], "\"")
			data[rsid] = geno
		}
	}
	return data, nil
}

// parseFTDNA2 parses the alternate FamilyTreeDNA CSV layout:
// rsid,chrom,pos,allele1,allele2.
func parseFTDNA2(filename string) (map[string]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() { _ = f.Close() }()

	data := make(map[string]string)
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 64<<10), 4<<20)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}

		parts := strings.Split(line, ",")
		if len(parts) < 5 {
			parts = strings.Fields(line)
		}

		if len(parts) == 5 {
			rsid := strings.Trim(parts[0], "\"")
			// allele1 + allele2
			g1 := strings.Trim(parts[3], "\"")
			g2 := strings.Trim(parts[4], "\"")
			data[rsid] = g1 + g2
		}
	}
	return data, nil
}
