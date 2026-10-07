package models

import (
	"strings"
	"testing"
)

func TestParseCoords(t *testing.T) {
	lines := []string{
		"rs1\t1\t1000\t1100\tA",
		"rs2\tX\t.\t2000\t.",
		"# comment",
		"notanrsid\t1\t1\t1\tA",
		"rs3\t2\t500\t.\tN", // invalid ref base
	}
	m := parseCoordsFromLines(lines...)
	if len(m) != 3 {
		t.Fatalf("len = %d, want 3", len(m))
	}
	c := m["rs1"]
	if c.Chrom != "1" || c.Pos19 != 1000 || c.Pos38 != 1100 || c.Ref != 'A' {
		t.Fatalf("rs1 = %+v", c)
	}
	c = m["rs2"]
	if c.Pos19 != 0 || c.Pos38 != 2000 || c.Ref != 0 {
		t.Fatalf("rs2 = %+v", c)
	}
	if m["rs3"].Ref != 0 {
		t.Fatalf("N ref should be 0")
	}
}

func TestBuildWantedAndExpand(t *testing.T) {
	m := &Model{
		SNPs: []string{"rs1", "rs2"},
		Coords: []*Coord{
			{RsID: "rs1", Chrom: "1", Pos38: 100, Pos19: 150, Ref: 'A'},
			{RsID: "rs2", Chrom: "2", Pos38: 500, Ref: 0},
		},
	}
	w := BuildWanted([]*Model{m})

	// rsID keys present
	if _, ok := w.Keys["rs1"]; !ok {
		t.Error("rs1 key missing")
	}
	// coordinate keys present with expected refs
	if ref := w.Keys["hg38:1:100"]; ref != 'A' {
		t.Errorf("hg38:1:100 ref = %q", ref)
	}
	if ref := w.Keys["hg19:1:150"]; ref != 'A' {
		t.Errorf("hg19:1:150 ref = %q", ref)
	}
	if _, ok := w.Keys["hg38:2:500"]; !ok {
		t.Error("hg38:2:500 key missing")
	}
	// hg19 absent for rs2
	if _, ok := w.Keys["hg19:2:500"]; ok {
		t.Error("unexpected hg19 key for rs2")
	}

	// Expand: range query on hg38 chr1
	got := w.Expand("hg38", "1", 50, 200)
	if len(got) != 1 || got[0].Key != "hg38:1:100" {
		t.Fatalf("expand = %+v", got)
	}
	if len(w.Expand("hg38", "1", 201, 999)) != 0 {
		t.Error("expected empty expand")
	}
	if len(w.Expand("hg19", "2", 1, 100000)) != 0 {
		t.Error("expected empty expand for other build")
	}
}

func TestWantedNilSafety(t *testing.T) {
	var w *Wanted
	if _, ok := w.Lookup("rs1"); ok {
		t.Error("nil Wanted must not report keys")
	}
	if len(w.Expand("hg38", "1", 0, 1000)) != 0 {
		t.Error("nil Wanted must expand to nothing")
	}
}

func TestCoordKey(t *testing.T) {
	c := &Coord{Chrom: "1", Pos19: 111, Pos38: 222}
	if got := c.Key("hg38"); got != "hg38:1:222" {
		t.Fatalf("hg38 key = %q", got)
	}
	if got := c.Key("hg19"); got != "hg19:1:111" {
		t.Fatalf("hg19 key = %q", got)
	}
}

func TestListModelsHasData(t *testing.T) {
	// The five models with embedded data must load.
	for _, name := range []string{"K12b", "K36", "K13M2", "HarappaWorld", "AncientNearEast13"} {
		m, err := LoadModel(name)
		if err != nil {
			t.Fatalf("LoadModel(%s): %v", name, err)
		}
		if len(m.SNPs) != len(m.Frequencies) {
			t.Errorf("%s: SNP/frequency count mismatch", name)
		}
		if len(m.Populations) == 0 {
			t.Errorf("%s: no populations", name)
		}
		if len(m.Coords) != len(m.SNPs) {
			t.Errorf("%s: coords length mismatch", name)
		}
	}
}

func TestLoadModelUnknown(t *testing.T) {
	if _, err := LoadModel("NoSuchModel"); err == nil {
		t.Fatal("expected error for unknown model")
	} else if !strings.Contains(err.Error(), "not found") {
		t.Fatalf("unexpected error: %v", err)
	}
}

func TestLoadAllRegisteredModels(t *testing.T) {
	// Every registered model must load with internally consistent data.
	for _, name := range ListModels() {
		m, err := LoadModel(name)
		if err != nil {
			t.Errorf("LoadModel(%s): %v", name, err)
			continue
		}
		if len(m.SNPs) != len(m.Frequencies) {
			t.Errorf("%s: %d SNPs but %d frequency rows", name, len(m.SNPs), len(m.Frequencies))
		}
		if len(m.Populations) == 0 {
			t.Errorf("%s: no populations", name)
		}
		if got := len(m.Frequencies[0]); got != len(m.Populations) {
			t.Errorf("%s: frequency row width %d, want %d (populations)", name, got, len(m.Populations))
		}
		if len(m.Coords) != len(m.SNPs) {
			t.Errorf("%s: coords length %d, want %d", name, len(m.Coords), len(m.SNPs))
		}
		for _, pop := range m.Populations {
			if pop.Name == "" {
				t.Errorf("%s: empty population name", name)
			}
		}
	}
}
