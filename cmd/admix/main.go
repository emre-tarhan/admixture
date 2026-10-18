package main

import (
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/emre-tarhan/admixture/pkg/calculator"
	"github.com/emre-tarhan/admixture/pkg/models"
	"github.com/emre-tarhan/admixture/pkg/parser"
	"github.com/emre-tarhan/admixture/pkg/render"
	"github.com/emre-tarhan/admixture/pkg/result"
)

func main() {
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "admix:", err)
		os.Exit(1)
	}
}

func run(args []string) error {
	fs := flag.NewFlagSet("admix", flag.ContinueOnError)
	filePtr := fs.String("f", "", "Raw genome data file (required)")
	modelsPtr := fs.String("m", "", "Admixture models, comma separated (e.g. K7b,K12b). Default: all")
	vendorPtr := fs.String("v", "23andme", "Vendor format (23andme, ancestry, ftdna, ftdna2, wegene, myheritage, livingdna; VCF/gVCF is auto-detected)")
	buildPtr := fs.String("b", "auto", "Reference build for VCF position matching (auto, hg19, hg38)")
	homRefPtr := fs.String("homref", "auto", "Homozygous-reference assumption for model SNPs missing from a variant-only VCF (auto, on, off)")
	outputPtr := fs.String("o", "", "Output file")
	jsonPtr := fs.Bool("json", false, "Output in JSON format")
	prettyPtr := fs.Bool("pretty", false, "Render results as bar charts in the terminal")
	tolerancePtr := fs.Float64("t", 1e-4, "Optimization tolerance")
	listModelsPtr := fs.Bool("models", false, "List available models and exit")

	fs.Usage = func() {
		out := fs.Output()
		fmt.Fprintf(out, `admix — admixture (ancestry composition) calculator

Usage:
  admix -f <raw-data-file> [flags]

Examples:
  admix -f genome.txt -m K12b -v 23andme     23andMe raw data, Dodecad K12b
  admix -f genome.vcf.gz -m K12b,K36         VCF (auto-detected), two models
  admix -f exome.vcf.gz -m K12b -homref off  exome VCF; do not assume hom-ref
  admix -f sample.g.vcf.gz -m K12b -pretty   gVCF, terminal bar charts
  admix -f genome.txt --json -o out.json     JSON results to a file
  admix -models                              list all available models

Flags:
`)
		fs.PrintDefaults()
	}

	if err := fs.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return nil // usage was printed by the flag package
		}
		return err
	}

	if *listModelsPtr {
		for _, name := range models.ListModels() {
			pops, err := models.GetPopulations(name)
			if err != nil {
				continue
			}
			fmt.Printf("%-18s %2d populations\n", name, len(pops))
		}
		return nil
	}

	if *jsonPtr && *prettyPtr {
		return fmt.Errorf("--json and -pretty are mutually exclusive")
	}
	homRef, err := parser.ParseRefMode(*homRefPtr)
	if err != nil {
		return err
	}

	// Validate args
	if *filePtr == "" {
		fs.Usage()
		return fmt.Errorf("no raw data file given (-f)")
	}

	// Determine models and load them up front: they define which rsIDs and
	// coordinates the VCF parser extracts.
	var modelList []string
	if *modelsPtr == "" {
		modelList = models.ListModels()
	} else {
		for _, m := range strings.Split(*modelsPtr, ",") {
			modelList = append(modelList, strings.TrimSpace(m))
		}
	}

	loaded := make([]*models.Model, 0, len(modelList))
	for _, modelName := range modelList {
		model, err := models.LoadModel(modelName)
		if err != nil {
			fmt.Printf("Error loading model %s: %v\n", modelName, err)
			continue
		}
		loaded = append(loaded, model)
	}
	if len(loaded) == 0 {
		return fmt.Errorf("no requested model could be loaded")
	}

	// Parse Raw Data (format is detected from file content; the -v flag
	// only matters for text files, so the resolved source is printed after)
	wanted := models.BuildWanted(loaded)
	fmt.Printf("Parsing raw data from %s...\n", *filePtr)
	rawData, src, err := parser.ParseFile(*filePtr, *vendorPtr, wanted, parser.Build(*buildPtr), homRef)
	if err != nil {
		return fmt.Errorf("parsing raw data: %w", err)
	}
	source := src.Vendor
	if source == "vcf" && src.Build != "" {
		source = "VCF, " + src.Build
	} else if source == "vcf" {
		source = "VCF, build undetermined"
	}
	fmt.Printf("Loaded %d genotypes (%s)", len(rawData), source)
	if src.Vendor == "vcf" {
		switch {
		case src.Build == "":
			fmt.Printf(" — rsID matching only; pass -b hg19 or -b hg38 for position matching")
		case src.Stats.HomRef == "assumed":
			fmt.Printf(" — rsID + position matching; %s unobserved model SNPs called homozygous-reference", comma(src.Stats.ImputedHomRef))
		case src.Stats.HomRef == "auto-skip":
			fmt.Printf(" — rsID + position matching; hom-ref not assumed for this variant-only file")
			fmt.Printf(" (smaller than a whole genome; pass -homref on to force)")
		default:
			fmt.Printf(" — rsID + position matching")
		}
	}
	fmt.Println()

	// Run Calculation for each model
	results := make([]result.Result, 0, len(loaded))
	for _, model := range loaded {
		fmt.Printf("Calculating %s...\n", model.Name)

		res, stats, err := calculator.Calculate(model, rawData, *tolerancePtr)
		if err != nil {
			fmt.Printf("Error calculating %s: %v\n", model.Name, err)
			continue
		}

		results = append(results, result.Result{
			Model:       model.Name,
			Proportions: res,
			MatchedSNPs: stats.MatchedSNPs,
			TotalSNPs:   stats.TotalSNPs,
		})
	}

	// Output
	if *jsonPtr {
		b, err := json.MarshalIndent(results, "", "  ")
		if err != nil {
			return err
		}
		return writeOutput(*outputPtr, b)
	}

	var out strings.Builder
	if *prettyPtr {
		render.Pretty(&out, results)
	} else {
		for _, r := range results {
			out.WriteString(fmt.Sprintf("Model: %s\n", r.Model))
			for _, e := range r.Sorted(0.0001) {
				out.WriteString(fmt.Sprintf("%s: %.2f%%\n", e.Name, e.Share*100))
			}
			out.WriteString(fmt.Sprintf("(%s of %s SNPs matched)\n\n",
				comma(r.MatchedSNPs), comma(r.TotalSNPs)))
		}
	}
	return writeOutput(*outputPtr, []byte(out.String()))
}

func writeOutput(outputPath string, data []byte) error {
	if outputPath == "" {
		_, _ = os.Stdout.Write(data)
		return nil
	}
	if err := os.WriteFile(outputPath, data, 0644); err != nil {
		return err
	}
	fmt.Printf("Results written to %s\n", outputPath)
	return nil
}

// comma formats n with thousands separators for progress lines.
func comma(n int) string {
	s := fmt.Sprint(n)
	if len(s) <= 3 {
		return s
	}
	var parts []string
	for len(s) > 3 {
		parts = append([]string{s[len(s)-3:]}, parts...)
		s = s[:len(s)-3]
	}
	return s + "," + strings.Join(parts, ",")
}
