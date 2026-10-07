package calculator

import (
	"github.com/emre-tarhan/admixture/pkg/models"
	"math"
	"sort"
)

// Stats carries diagnostics about one calculation run.
type Stats struct {
	MatchedSNPs int // model SNPs found in the raw data
	TotalSNPs   int // model SNPs overall
}

// Calculate performs the admixture analysis.
func Calculate(model *models.Model, rawData map[string]string, tolerance float64) (map[string]float64, Stats, error) {
	// 1. Genotype Matches
	gMajor, gMinor, stats := genotypeMatches(model, rawData)

	// 2. Optimization
	nPop := len(model.Populations)

	// Initial guess: uniform
	q := make([]float64, nPop)
	for i := range q {
		q[i] = 1.0 / float64(nPop)
	}

	// Optimization Loop (Projected Gradient Descent with Adaptive LR)
	maxIter := 5000
	learningRate := 0.001 // Start with a reasonable rate

	currentLoss, grad := calculateLossAndGradient(q, model.Frequencies, gMajor, gMinor)

	for iter := 0; iter < maxIter; iter++ {
		// Propose new Q
		qNew := make([]float64, nPop)
		for i := 0; i < nPop; i++ {
			qNew[i] = q[i] - learningRate*grad[i]
		}
		qNew = projectToSimplex(qNew)

		// Calculate loss at new Q
		lossNew, gradNew := calculateLossAndGradient(qNew, model.Frequencies, gMajor, gMinor)

		if lossNew < currentLoss {
			// Accept step
			if math.Abs(currentLoss-lossNew) < tolerance {
				q = qNew
				break
			}
			q = qNew
			currentLoss = lossNew
			grad = gradNew
			learningRate *= 1.1 // Boldly go
		} else {
			// Reject step
			learningRate *= 0.5
			if learningRate < 1e-10 {
				break
			}
		}
	}

	// Format results
	results := make(map[string]float64)
	for i, pop := range model.Populations {
		results[pop.Name] = q[i]
	}

	return results, stats, nil
}

// genotypeMatches counts major/minor allele copies per model SNP, looking
// each SNP up by rsID first and then by position (either build).
func genotypeMatches(model *models.Model, rawData map[string]string) ([]float64, []float64, Stats) {
	nSNPs := len(model.SNPs)
	gMajor := make([]float64, nSNPs)
	gMinor := make([]float64, nSNPs)
	stats := Stats{TotalSNPs: nSNPs}

	for i, snp := range model.SNPs {
		genotype, ok := rawData[snp]
		if !ok && i < len(model.Coords) && model.Coords[i] != nil {
			// Fall back to positional keys (VCF/gVCF input without rsIDs).
			c := model.Coords[i]
			if c.Pos38 > 0 {
				genotype, ok = rawData[c.Key("hg38")]
			}
			if !ok && c.Pos19 > 0 {
				genotype, ok = rawData[c.Key("hg19")]
			}
		}
		if !ok {
			// Missing data: contributes to neither loss nor gradient
			// (same as vendor files with holes).
			continue
		}
		stats.MatchedSNPs++

		// Check first allele
		a1 := string(genotype[0])
		// Check second allele (if exists)
		a2 := ""
		if len(genotype) > 1 {
			a2 = string(genotype[len(genotype)-1])
		} else {
			a2 = a1 // Hemizygous: a single allele counts as both copies.
		}

		major := model.MajorAlleles[i]
		minor := model.MinorAlleles[i]

		if a1 == major {
			gMajor[i] += 1
		}
		if a2 == major {
			gMajor[i] += 1
		}

		if a1 == minor {
			gMinor[i] += 1
		}
		if a2 == minor {
			gMinor[i] += 1
		}
	}
	return gMajor, gMinor, stats
}

func calculateLossAndGradient(q []float64, freqs [][]float64, gMajor, gMinor []float64) (float64, []float64) {
	loss := 0.0
	grad := make([]float64, len(q))
	nPop := len(q)

	// Parallelize this loop if performance is an issue.
	for i := range gMajor {
		// Skip if no data for this SNP
		if gMajor[i] == 0 && gMinor[i] == 0 {
			continue
		}

		// P_major = sum(F_ik * q_k)
		pMajor := 0.0
		for k := 0; k < nPop; k++ {
			pMajor += freqs[i][k] * q[k]
		}

		// Numerical stability
		if pMajor < 1e-9 {
			pMajor = 1e-9
		}
		if pMajor > 1.0-1e-9 {
			pMajor = 1.0 - 1e-9
		}

		pMinor := 1.0 - pMajor

		// Loss = - (gMajor * log(pMajor) + gMinor * log(pMinor))
		loss -= (gMajor[i]*math.Log(pMajor) + gMinor[i]*math.Log(pMinor))

		// Gradient
		// d(-L)/dq_k = - [ gMajor * F_ik / pMajor - gMinor * F_ik / pMinor ]
		//            = - [ gMajor/pMajor - gMinor/pMinor ] * F_ik

		term := (gMajor[i] / pMajor) - (gMinor[i] / pMinor)
		for k := 0; k < nPop; k++ {
			grad[k] -= term * freqs[i][k]
		}
	}
	return loss, grad
}

// projectToSimplex projects vector v onto the simplex constraints (sum=1, non-negative).
// Reference: https://eng.ucmerced.edu/people/wwang5/papers/SimplexProj.pdf
func projectToSimplex(v []float64) []float64 {
	n := len(v)
	u := make([]float64, n)
	copy(u, v)

	// Sort u in descending order
	sort.Sort(sort.Reverse(sort.Float64Slice(u)))

	rho := 0
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += u[i]
		if u[i]+(1.0-sum)/float64(i+1) > 0 {
			rho = i
		}
	}

	sumRho := 0.0
	for i := 0; i <= rho; i++ {
		sumRho += u[i]
	}
	theta := (1.0 - sumRho) / float64(rho+1)

	w := make([]float64, n)
	for i := 0; i < n; i++ {
		val := v[i] + theta
		if val < 0 {
			val = 0
		}
		w[i] = val
	}
	return w
}
