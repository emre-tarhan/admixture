package calculator

import (
	"math"
	"testing"

	"github.com/emre-tarhan/admixture/pkg/models"
)

// twoPopModel builds a tiny synthetic model. NOTE on conventions: the .F
// frequency matrix holds the frequency of the MAJOR allele (verified
// empirically — values routinely exceed 0.5, which is impossible for minor
// frequencies). Population A has major-allele frequency 0.3, B 0.9.
func twoPopModel(n int) *models.Model {
	m := &models.Model{
		Name:        "test",
		Populations: []models.Population{{Name: "A"}, {Name: "B"}},
	}
	for i := 0; i < n; i++ {
		m.SNPs = append(m.SNPs, "rs"+string(rune('a'+i)))
		m.MinorAlleles = append(m.MinorAlleles, "T")
		m.MajorAlleles = append(m.MajorAlleles, "C")
		m.Frequencies = append(m.Frequencies, []float64{0.4, 0.8})
	}
	return m
}

func TestCalculateRecoversMixture(t *testing.T) {
	m := twoPopModel(500)

	// An individual with major-allele (C) rate 0.5: TT at 20% of sites,
	// TC at 60%, CC at 20% → (2·0 + 6·1 + 2·2) = 10 major copies over
	// 20 alleles = 0.5.
	raw := map[string]string{}
	for i, snp := range m.SNPs {
		switch {
		case i%10 < 2:
			raw[snp] = "TT"
		case i%10 < 8:
			raw[snp] = "TC"
		default:
			raw[snp] = "CC"
		}
	}

	res, stats, err := Calculate(m, raw, 1e-4)
	if err != nil {
		t.Fatal(err)
	}
	if stats.MatchedSNPs != 500 || stats.TotalSNPs != 500 {
		t.Fatalf("stats = %+v", stats)
	}
	sum := res["A"] + res["B"]
	if math.Abs(sum-1) > 1e-6 {
		t.Fatalf("proportions sum to %f", sum)
	}
	// MLE: p̂ = 0.5 = 0.4·A + 0.8·(1−A) → A = 0.75.
	if res["A"] < 0.6 || res["A"] > 0.9 {
		t.Fatalf("expected A ≈ 0.75, got %v", res)
	}
}

func TestCalculateMissingSNPsCounted(t *testing.T) {
	m := twoPopModel(100)
	raw := map[string]string{"rsa": "TT", "rsb": "CC"} // only 2 of 100
	_, stats, err := Calculate(m, raw, 1e-4)
	if err != nil {
		t.Fatal(err)
	}
	if stats.MatchedSNPs != 2 || stats.TotalSNPs != 100 {
		t.Fatalf("stats = %+v", stats)
	}
}

func TestCalculatePositionalKeys(t *testing.T) {
	// VCF-style rawData: lookup by coordinate key attached to the model.
	m := twoPopModel(3)
	m.Coords = []*models.Coord{
		{RsID: "rsa", Chrom: "1", Pos38: 100},
		{RsID: "rsb", Chrom: "1", Pos38: 200},
		{RsID: "rsc", Chrom: "1", Pos38: 300},
	}
	raw := map[string]string{
		"hg38:1:100": "TT",
		"hg38:1:200": "TT",
		"hg38:1:300": "CC",
	}
	_, stats, err := Calculate(m, raw, 1e-4)
	if err != nil {
		t.Fatal(err)
	}
	if stats.MatchedSNPs != 3 {
		t.Fatalf("positional matching failed: %+v", stats)
	}
}

func TestProjectToSimplex(t *testing.T) {
	cases := []struct {
		in   []float64
		want []float64
	}{
		{[]float64{3, 1, 2}, []float64{1, 0, 0}},
		{[]float64{0.5, 0.5}, []float64{0.5, 0.5}},
		{[]float64{-1, 2}, []float64{0, 1}},
	}
	for i, c := range cases {
		got := projectToSimplex(c.in)
		sum := 0.0
		for _, v := range got {
			if v < 0 {
				t.Errorf("case %d: negative component %f", i, v)
			}
			sum += v
		}
		if math.Abs(sum-1) > 1e-9 {
			t.Errorf("case %d: sum = %f", i, sum)
		}
		for j := range c.want {
			if math.Abs(got[j]-c.want[j]) > 1e-9 {
				t.Errorf("case %d: got %v, want %v", i, got, c.want)
				break
			}
		}
	}
}

func TestLossGradientFinite(t *testing.T) {
	m := twoPopModel(10)
	q := []float64{0.6, 0.4}
	gMajor := []float64{1, 2, 1, 0, 2, 1, 1, 2, 0, 1}
	gMinor := []float64{1, 0, 1, 2, 0, 1, 1, 0, 2, 1}
	loss, grad := calculateLossAndGradient(q, m.Frequencies, gMajor, gMinor)
	if math.IsNaN(loss) || math.IsInf(loss, 0) {
		t.Fatalf("loss not finite: %f", loss)
	}
	for k, g := range grad {
		if math.IsNaN(g) || math.IsInf(g, 0) {
			t.Fatalf("grad[%d] not finite: %f", k, g)
		}
	}
}
