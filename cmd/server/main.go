package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"strconv"
	"strings"

	"github.com/emre-tarhan/admixture/pkg/calculator"
	"github.com/emre-tarhan/admixture/pkg/models"
	"github.com/emre-tarhan/admixture/pkg/parser"
	"github.com/emre-tarhan/admixture/pkg/result"
)

// DefaultPort is the port the API server listens on when neither -port nor
// PORT is given. Chosen to be distinctive (keys "ADMI" on a phone keypad)
// so it is unlikely to collide with anything else running locally.
const DefaultPort = 2364

func main() {
	portFlag := flag.Int("port", 0, "Port to listen on (overrides the PORT environment variable)")
	addrFlag := flag.String("addr", "127.0.0.1", "Address to bind (use 0.0.0.0 to expose — only behind an authenticating proxy)")
	flag.Parse()

	port := DefaultPort
	if env := os.Getenv("PORT"); env != "" {
		if p, err := strconv.Atoi(env); err == nil {
			port = p
		}
	}
	if *portFlag != 0 {
		port = *portFlag
	}

	http.HandleFunc("/health", handleHealth)
	http.HandleFunc("/calculate", handleCalculate)
	http.HandleFunc("/models", handleListModels)

	listen := *addrFlag + ":" + strconv.Itoa(port)
	fmt.Printf("Server listening on %s...\n", listen)
	log.Fatal(http.ListenAndServe(listen, nil))
}

func handleHealth(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]string{"status": "ok"})
}

func handleListModels(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]interface{}{
		"models": models.ListModels(),
	})
}

func handleCalculate(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	// Parse multipart form
	err := r.ParseMultipartForm(10 << 20) // 10 MB limit
	if err != nil {
		http.Error(w, "Unable to parse form", http.StatusBadRequest)
		return
	}

	file, _, err := r.FormFile("raw_data")
	if err != nil {
		http.Error(w, "Missing raw_data file", http.StatusBadRequest)
		return
	}
	defer func() { _ = file.Close() }()

	vendor := r.FormValue("vendor")
	if vendor == "" {
		vendor = "23andme" // Default
	}

	modelName := r.FormValue("model")
	if modelName == "" {
		modelName = "K7b" // Default
	}

	build := r.FormValue("build")
	if build == "" {
		build = "auto"
	}

	// Save to temp file
	tempFile, err := os.CreateTemp("", "admix-upload-*.txt")
	if err != nil {
		http.Error(w, "Server error", http.StatusInternalServerError)
		return
	}
	defer func() {
		_ = tempFile.Close() // best-effort cleanup of a temp file
		_ = os.Remove(tempFile.Name())
	}()

	if _, err = io.Copy(tempFile, file); err != nil {
		http.Error(w, "Server error", http.StatusInternalServerError)
		return
	}

	// Load models first: they define which rsIDs and coordinates the VCF
	// parser extracts.
	modelNames := strings.Split(modelName, ",")
	loaded := make([]*models.Model, 0, len(modelNames))
	for _, m := range modelNames {
		m = strings.TrimSpace(m)
		model, err := models.LoadModel(m)
		if err != nil {
			continue // unknown or undataed models are skipped silently
		}
		loaded = append(loaded, model)
	}
	if len(loaded) == 0 {
		http.Error(w, "No requested model is available", http.StatusBadRequest)
		return
	}

	if _, err = tempFile.Seek(0, io.SeekStart); err != nil {
		http.Error(w, "Server error", http.StatusInternalServerError)
		return
	}
	rawData, _, err := parser.ParseFile(tempFile.Name(), vendor, models.BuildWanted(loaded), parser.Build(build))
	if err != nil {
		http.Error(w, fmt.Sprintf("Error parsing data: %v", err), http.StatusBadRequest)
		return
	}

	results := make([]result.Result, 0, len(loaded))
	for _, model := range loaded {
		res, stats, err := calculator.Calculate(model, rawData, 1e-4)
		if err != nil {
			continue
		}
		results = append(results, result.Result{
			Model:       model.Name,
			Proportions: res,
			MatchedSNPs: stats.MatchedSNPs,
			TotalSNPs:   stats.TotalSNPs,
		})
	}

	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(results)
}
