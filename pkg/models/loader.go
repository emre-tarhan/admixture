package models

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"strconv"
	"strings"
)

// LoadModel loads the model data from the embedded filesystem.
func LoadModel(name string) (*Model, error) {
	populations, err := GetPopulations(name)
	if err != nil {
		return nil, err
	}
	nPop := len(populations)

	// Load SNPs
	// Path relative to pkg/models
	snpFile := "data/" + name + ".alleles"
	f, err := content.Open(snpFile)
	if err != nil {
		return nil, fmt.Errorf("failed to open allele file %s: %w", snpFile, err)
	}
	defer func() { _ = f.Close() }()

	var snps, minor, major []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		parts := strings.Fields(line)
		if len(parts) < 3 {
			continue
		}
		snps = append(snps, parts[0])
		minor = append(minor, parts[1])
		major = append(major, parts[2])
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading allele file: %w", err)
	}

	// Frequency matrices are stored gzip-compressed (<Model>.<N>.F.gz) to
	// keep the binary small; they decompress transparently on load.
	freqFile := fmt.Sprintf("data/%s.%d.F.gz", name, nPop)
	f2, err := content.Open(freqFile)
	if err != nil {
		return nil, fmt.Errorf("failed to open frequency file %s: %w", freqFile, err)
	}
	defer func() { _ = f2.Close() }()
	gz, err := gzip.NewReader(f2)
	if err != nil {
		return nil, fmt.Errorf("failed to decompress frequency file %s: %w", freqFile, err)
	}
	defer func() { _ = gz.Close() }()

	var freqs [][]float64
	scanner2 := bufio.NewScanner(gz)
	for scanner2.Scan() {
		line := strings.TrimSpace(scanner2.Text())
		if line == "" {
			continue
		}
		parts := strings.Fields(line)
		var row []float64
		for _, val := range parts {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid frequency value '%s': %w", val, err)
			}
			row = append(row, f)
		}
		if len(row) != nPop {
			return nil, fmt.Errorf("frequency row has %d columns, expected %d", len(row), nPop)
		}
		freqs = append(freqs, row)
	}
	if err := scanner2.Err(); err != nil {
		return nil, fmt.Errorf("error reading frequency file: %w", err)
	}

	// Validation
	if len(snps) != len(freqs) {
		return nil, fmt.Errorf("mismatch between SNP count (%d) and Frequency row count (%d)", len(snps), len(freqs))
	}

	// Attach reference-build coordinates where the embedded index has them
	// (used for position-based VCF/gVCF matching; rsID-only when absent).
	coords := make([]*Coord, len(snps))
	for i, snp := range snps {
		coords[i] = coordFor(snp)
	}

	return &Model{
		Name:         name,
		Populations:  populations,
		SNPs:         snps,
		MinorAlleles: minor,
		MajorAlleles: major,
		Frequencies:  freqs,
		Coords:       coords,
	}, nil
}
