package parser

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/emre-tarhan/admixture/pkg/models"
)

// Build is a human genome reference build.
type Build string

const (
	BuildAuto Build = "auto"
	Build19   Build = "hg19"
	Build38   Build = "hg38"
)

// VCFStats carries diagnostics from a VCF/gVCF parse run.
type VCFStats struct {
	Records       int // data lines seen
	Kept          int // records that produced at least one genotype key
	Filtered      int // records skipped due to a failing FILTER column
	BlockCoords   int // genotype keys emitted from gVCF reference blocks
	ImputedHomRef int // model SNPs called homozygous-reference because the VCF had no record
	// HomRef reports what happened to model SNPs the file does not list:
	// "off", "gvcf" (reference blocks already encode hom-ref), "auto-skip"
	// (not a whole-genome variant VCF), or "assumed".
	HomRef      string
	BuildSource string
}

// A variant-only whole-genome VCF lists several million sites spread across
// the autosomes. An exome is one to two orders of magnitude smaller, and a
// single chromosome slice does not cover 20 autosomes. Both of those must
// keep "no record" meaning "not genotyped" — filling them as hom-ref would
// invent calls outside the assayed regions.
const (
	wgsMinRecords   = 1_000_000
	wgsMinAutosomes = 20
)

// vcfBuildSignatures maps known chromosome lengths (1-based contig
// declarations) to reference builds. A few chromosomes per build are enough
// — they never overlap between hg19 and hg38.
var vcfBuildSignatures = map[int]Build{
	249250621: Build19, // chr1
	243199373: Build19, // chr2
	198022430: Build19, // chr3
	155270560: Build19, // chrX
	248956422: Build38, // chr1
	242193529: Build38, // chr2
	198295559: Build38, // chr3
	156040895: Build38, // chrX
}

// ParseVCF parses a (possibly bgzip/gzip-compressed) VCF or gVCF file into
// lookup keys → genotype, where keys are rsIDs and/or
// "hgNN:chrom:pos" coordinates drawn from the wanted set.
//
// Matching strategy:
//   - Records with a dbSNP ID emit their rsID key (works for any build).
//   - When the reference build is known, records additionally emit their
//     coordinate key (covers WGS callsets that carry no rsIDs at all).
//   - gVCF reference blocks (ALT = <*>/<NON_REF>, INFO END=) assert the
//     sample is homozygous-reference across the span: for every wanted
//     coordinate inside the span with a known expected reference allele,
//     the reference base is emitted as a hom-ref genotype.
//
// wanted may be nil, in which case every rsID is kept (coordinate keys are
// still only emitted when the build is known). Indels, MNPs, symbolic
// alleles, non-PASS records and no-calls (GT ./. ) are skipped.
//
// On a variant-only whole-genome VCF (ref mode "auto", or "on" for any
// variant-only file) a model SNP on a primary chromosome that produced at
// least one wanted hit, with no record of its own, is emitted as
// homozygous for the index reference allele. gVCF reference blocks already
// carry those calls, so a file that contains any block is left alone.
// Sites that do have a record — including filtered, indel and no-call
// records — are never replaced.
//
// Returns the parsed genotypes, the build used for coordinate matching
// ("" when undetermined) and per-run statistics.
func ParseVCF(filename string, wanted *models.Wanted, force Build, ref RefMode) (map[string]string, Build, VCFStats, error) {
	var stats VCFStats

	f, err := os.Open(filename)
	if err != nil {
		return nil, "", stats, err
	}
	defer func() { _ = f.Close() }()

	raw, err := decompressMaybe(f)
	if err != nil {
		return nil, "", stats, err
	}
	defer func() { _ = raw.Close() }()

	scanner := bufio.NewScanner(raw)
	scanner.Buffer(make([]byte, 0, 256<<10), 8<<20)

	// ---- header: contig lengths for build detection, column layout ----
	sampleCol := -1
	detected := Build("")
headerLoop:
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "##contig=<") {
			if detected == "" && force == BuildAuto {
				detected = detectBuildFromContig(line)
			}
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			cols := strings.Split(line, "\t")
			for i, c := range cols {
				if c == "FORMAT" {
					sampleCol = i + 1
				}
			}
			break headerLoop
		}
	}
	if sampleCol < 0 {
		return nil, "", stats, fmt.Errorf("no sample genotype columns found (sites-only VCF?)")
	}
	build := force
	buildSource := "forced"
	if build == BuildAuto {
		build, buildSource = detected, "contig lengths"
	}
	if build == BuildAuto || build == "" {
		build = ""
		buildSource = "undetermined"
	}
	stats.BuildSource = buildSource

	data := make(map[string]string, 256<<10)
	gtIndex, lastFormat := -1, ""
	// observed marks wanted coordinates that have a record, usable or not,
	// so imputation cannot overwrite a filtered or non-SNP call.
	observed := make(map[string]struct{})
	wantedHits := make(map[string]int)
	perChrom := make(map[string]int)
	sawRefBlock := false

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}
		stats.Records++
		fields := strings.Split(line, "\t")
		if len(fields) <= sampleCol {
			continue
		}

		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}
		chrom := normalizeChrom(fields[0])
		if isPrimaryChrom(chrom) {
			perChrom[chrom]++
		}
		if build != "" && wanted != nil {
			key := coordKey(string(build), chrom, pos)
			if _, ok := wanted.Lookup(key); ok {
				observed[key] = struct{}{}
				wantedHits[chrom]++
			}
		}

		// FILTER: keep "." (unfiltered) and PASS. The site stays observed,
		// so it is not later filled in as homozygous-reference.
		if filt := fields[6]; filt != "." && filt != "PASS" {
			stats.Filtered++
			continue
		}

		ref := fields[3]
		id := fields[2]
		format := fields[8]
		sample := fields[sampleCol]
		// Any reference-block record marks the file as a gVCF, even when
		// this particular block is a no-call and emits nothing.
		if isRefBlock(fields[4], ref) {
			sawRefBlock = true
		}

		if format != lastFormat {
			gtIndex = -1
			for i, key := range strings.Split(format, ":") {
				if key == "GT" {
					gtIndex = i
					break
				}
			}
			lastFormat = format
		}
		if gtIndex < 0 {
			continue // no genotype field in this record
		}
		sampleFields := strings.Split(sample, ":")
		if gtIndex >= len(sampleFields) {
			continue
		}
		gt := sampleFields[gtIndex]

		if isRefBlock(fields[4], ref) {
			// gVCF reference block: sample is hom-ref across [pos, END]
			// when the GT is fully called 0/0.
			if !isHomRef(gt) {
				continue // no-call or partial: nothing to emit
			}
			end := blockEnd(fields[7], pos)
			n := expandRefBlock(data, wanted, string(build), chrom, pos, end)
			stats.BlockCoords += n
			if n > 0 {
				stats.Kept++
			}
			continue
		}

		genotype, ok := resolveGenotype(gt, ref, fields[4])
		if !ok {
			continue
		}

		kept := false
		if id != "." {
			for _, rsid := range strings.Split(id, ";") {
				if !strings.HasPrefix(rsid, "rs") {
					continue
				}
				if _, exists := data[rsid]; !exists {
					data[rsid] = genotype
					kept = true
				}
			}
		}
		if build != "" {
			key := coordKey(string(build), chrom, pos)
			if _, exists := data[key]; !exists {
				data[key] = genotype
				kept = true
			}
		}
		if kept {
			stats.Kept++
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, build, stats, err
	}
	stats.ImputedHomRef, stats.HomRef = imputeMissingHomRef(data, wanted, observed, wantedHits, perChrom, build, ref, sawRefBlock, stats.Records)
	return data, build, stats, nil
}

// imputeMissingHomRef fills wanted coordinates the VCF never mentioned.
// See ParseVCF for when that is the right reading of a variant-only file.
func imputeMissingHomRef(data map[string]string, wanted *models.Wanted, observed map[string]struct{}, wantedHits, perChrom map[string]int, build Build, ref RefMode, sawRefBlock bool, records int) (int, string) {
	if ref == "" {
		ref = RefAuto
	}
	if ref == RefOff || wanted == nil || build == "" {
		return 0, "off"
	}
	if sawRefBlock {
		return 0, "gvcf"
	}
	if ref == RefAuto && !looksLikeWGS(records, perChrom) {
		return 0, "auto-skip"
	}
	prefix := string(build) + ":"
	n := 0
	for key, refBase := range wanted.Keys {
		if refBase == 0 || !strings.HasPrefix(key, prefix) {
			continue
		}
		if _, exists := data[key]; exists {
			continue
		}
		if _, seen := observed[key]; seen {
			continue
		}
		chrom := chromOfKey(key)
		// A chromosome with no wanted hit at all is either absent from the
		// file or on the wrong reference build. Filling it would turn every
		// model SNP into a reference homozygote.
		if !isPrimaryChrom(chrom) || wantedHits[chrom] == 0 {
			continue
		}
		data[key] = string([]byte{refBase, refBase})
		n++
	}
	return n, "assumed"
}

// looksLikeWGS reports whether a variant-only VCF has the size and
// chromosomal spread of a whole-genome callset.
func looksLikeWGS(records int, perChrom map[string]int) bool {
	if records < wgsMinRecords {
		return false
	}
	autosomes := 0
	for i := 1; i <= 22; i++ {
		if perChrom[strconv.Itoa(i)] > 0 {
			autosomes++
		}
	}
	return autosomes >= wgsMinAutosomes
}

func isPrimaryChrom(chrom string) bool {
	switch chrom {
	case "X", "Y", "M":
		return true
	}
	n, err := strconv.Atoi(chrom)
	return err == nil && n >= 1 && n <= 22
}

func chromOfKey(key string) string {
	parts := strings.Split(key, ":")
	if len(parts) < 3 {
		return ""
	}
	return parts[1]
}

func expandRefBlock(data map[string]string, wanted *models.Wanted, build, chrom string, start, end int) int {
	coords := wanted.Expand(build, chrom, start, end)
	emitted := 0
	for _, c := range coords {
		if c.Ref == 0 {
			continue
		}
		if _, exists := data[c.Key]; !exists {
			data[c.Key] = string(c.Ref) + string(c.Ref)
			emitted++
		}
	}
	return emitted
}

// isRefBlock reports whether a record is a gVCF reference-block record:
// no real ALT allele (missing, <*> or <NON_REF> only).
func isRefBlock(altField, ref string) bool {
	if altField == "." || altField == "" {
		return true
	}
	for _, alt := range strings.Split(altField, ",") {
		if alt != "" && !strings.HasPrefix(alt, "<") && alt != "*" {
			return false
		}
	}
	return true
}

// isHomRef reports whether a GT is fully called and contains only allele 0
// (e.g. "0/0", "0|0", "0").
func isHomRef(gt string) bool {
	tokens := splitGT(gt)
	if len(tokens) == 0 {
		return false
	}
	for _, token := range tokens {
		if token != "0" {
			return false
		}
	}
	return true
}

// blockEnd extracts END= from INFO, defaulting to pos.
func blockEnd(info string, pos int) int {
	const tag = "END="
	for _, kv := range strings.Split(info, ";") {
		if strings.HasPrefix(kv, tag) {
			if n, err := strconv.Atoi(kv[len(tag):]); err == nil {
				return n
			}
		}
	}
	return pos
}

// resolveGenotype converts a GT (e.g. "0|1", "1/1", "2") into a genotype
// string of allele bases ("AT", "TT", "T"). It returns false when the
// record is a no-call, partially called, or involves an allele that is not
// a single base (indels, MNPs, symbolic and spanning-deletion alleles).
func resolveGenotype(gt, ref, altField string) (string, bool) {
	alts := strings.Split(altField, ",")
	altBase := func(idx int) (byte, bool) {
		if idx < 0 || idx >= len(alts) {
			return 0, false
		}
		a := alts[idx]
		if len(a) == 1 && isVCFBase(a[0]) {
			return a[0], true
		}
		return 0, false
	}
	if len(ref) != 1 || !isVCFBase(ref[0]) {
		return "", false
	}
	var b strings.Builder
	for _, token := range splitGT(gt) {
		if token == "." || token == "" {
			return "", false
		}
		n, err := strconv.Atoi(token)
		if err != nil {
			return "", false
		}
		if n == 0 {
			b.WriteByte(ref[0])
			continue
		}
		allele, ok := altBase(n - 1)
		if !ok {
			return "", false
		}
		b.WriteByte(allele)
	}
	return b.String(), b.Len() > 0
}

func splitGT(gt string) []string {
	// Keep empty tokens: "./." and "0/." are partial no-calls that must
	// not silently become valid haploid genotypes.
	return strings.Split(strings.ReplaceAll(gt, "|", "/"), "/")
}

func isVCFBase(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		return true
	}
	return false
}

// normalizeChrom maps VCF chromosome names to the canonical form used by
// the coordinate index: no "chr" prefix, 23 -> X, 24 -> Y, MT -> M.
func normalizeChrom(chrom string) string {
	c := strings.TrimPrefix(chrom, "chr")
	switch c {
	case "23":
		return "X"
	case "24":
		return "Y"
	case "MT", "mt":
		return "M"
	}
	return c
}

func coordKey(build string, chrom string, pos int) string {
	return build + ":" + chrom + ":" + strconv.Itoa(pos)
}

func decompressMaybe(f *os.File) (io.ReadCloser, error) {
	br := bufio.NewReader(f)
	magic, err := br.Peek(2)
	if err != nil {
		return nil, err
	}
	if magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(br)
		if err != nil {
			return nil, err
		}
		return gz, nil
	}
	// wrap so Close stays consistent and the buffered head is not lost
	return &nopCloser{Reader: br}, nil
}

type nopCloser struct{ *bufio.Reader }

func (n *nopCloser) Close() error { return nil }

// detectBuildFromContig looks for a known chromosome length inside a
// "##contig=<ID=chr1,length=...>" header line.
func detectBuildFromContig(line string) Build {
	for _, kv := range strings.Split(strings.Trim(strings.TrimPrefix(line, "##contig=<"), ">"), ",") {
		if v, ok := strings.CutPrefix(kv, "length="); ok {
			length, _ := strconv.Atoi(v)
			return vcfBuildSignatures[length]
		}
	}
	return ""
}
