package models

// Population represents a population in an admixture model.
type Population struct {
	Name   string // English name
	NameZH string // localized display name (Chinese / Turkish)
}

// Model represents a loaded admixture model.
type Model struct {
	Name         string
	Populations  []Population
	SNPs         []string
	MinorAlleles []string
	MajorAlleles []string
	Frequencies  [][]float64 // [SNP Index][Population Index]
	// Coords[i] holds reference-build coordinates for SNPs[i], or nil when
	// the SNP has no entry in the embedded coordinate index. VCF/gVCF input
	// is matched positionally through these (see pkg/parser).
	Coords []*Coord
}
