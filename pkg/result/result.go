// Package result defines the shared admixture result shape used by the
// CLI, the HTTP API and the terminal renderer.
package result

// Result holds the outcome of one model's calculation.
type Result struct {
	Model       string             `json:"model"`
	Proportions map[string]float64 `json:"proportions"` // population -> share of ancestry (0-1)
	MatchedSNPs int                `json:"matched_snps"`
	TotalSNPs   int                `json:"total_snps"`
}

// Sorted returns the populations ordered by proportion, descending.
// Contributions below minShare are dropped (0 keeps everything).
func (r Result) Sorted(minShare float64) []Entry {
	entries := make([]Entry, 0, len(r.Proportions))
	for name, share := range r.Proportions {
		if share < minShare {
			continue
		}
		entries = append(entries, Entry{Name: name, Share: share})
	}
	sortByShare(entries)
	return entries
}

// Entry is one population's share inside a Result.
type Entry struct {
	Name  string
	Share float64
}

func sortByShare(entries []Entry) {
	for i := 1; i < len(entries); i++ {
		for j := i; j > 0 && entries[j].Share > entries[j-1].Share; j-- {
			entries[j], entries[j-1] = entries[j-1], entries[j]
		}
	}
}
