# admixture

A fast, dependency-free admixture calculator in Go. It takes raw autosomal
DNA data — vendor text exports **or VCF/gVCF files** — and estimates your
ancestral composition as a mixture of reference populations, using
well-known community calculators such as **Dodecad K12b**, **Eurogenes
K36**, **HarappaWorld**, **AncientNearEast13** and **K13M2**.

Your genome file never leaves your machine: all model data is embedded in a
single static binary, and the optional API server binds to `127.0.0.1` only.

[![CI](https://github.com/emre-tarhan/admixture/actions/workflows/ci.yml/badge.svg)](https://github.com/emre-tarhan/admixture/actions/workflows/ci.yml)
[![Release](https://github.com/emre-tarhan/admixture/actions/workflows/release.yml/badge.svg)](https://github.com/emre-tarhan/admixture/actions/workflows/release.yml)

## How it works

1. **Parse** — reads your raw data file (23andMe-style text, VCF or gVCF)
   into genotypes keyed by rsID and/or chromosomal position.
2. **Match** — for every SNP in the selected model, counts how many of your
   two alleles match the model's major and minor allele, and reports how
   many of the model's SNPs were found in your input.
3. **Optimize** — finds the population mixture `q` (non-negative, summing to
   1) that maximizes the binomial likelihood of those counts, via
   projected gradient descent on the probability simplex.

Details: [docs/architecture.md](docs/architecture.md).

## Features

- **Single binary** — model files and the hg19/hg38 coordinate index are
  embedded with `go:embed`; nothing to download or configure.
- **Zero dependencies** — Go standard library only.
- **Fast** — a full K12b fit (166,770 SNPs) takes well under a second; a
  4M-record whole-genome VCF parses in ~20 s.
- **VCF/gVCF input** — works with or without rsID annotations, on hg19 and
  hg38 (auto-detected), including gVCF reference blocks. A variant-only
  whole-genome VCF is read as homozygous-reference at model SNPs it does
  not list.
- **CLI and HTTP API** — use it interactively or as a local service.
- **Multiple vendor formats** — 23andMe, AncestryDNA, FTDNA, MyHeritage,
  WeGene, LivingDNA.
- **Tested and linted** — race-detector test suite on Linux/Windows/macOS,
  golangci-lint, GoReleaser releases.

## Quick start

Download a binary from the [releases](https://github.com/emre-tarhan/admixture/releases),
or build from source (Go 1.25+):

```bash
go build -o admix ./cmd/admix
```

Run against a raw data file:

```bash
./admix -f genome.txt -m K12b -v 23andme          # vendor text export
./admix -f sample.vcf.gz -m K12b                  # VCF (auto-detected)
./admix -f sample.g.vcf.gz -m K12b,K36            # gVCF, two models
```

```
Parsing raw data from sample.vcf.gz...
Loaded 4177839 genotypes (VCF, hg38) — rsID + position matching; 95,282 unobserved model SNPs called homozygous-reference
Calculating K12b...
Model: K12b
Caucasus: 43.01%
Gedrosia: 15.65%
Southwest Asian: 12.45%
...
(164,974 of 166,770 SNPs matched)
```

Add `-pretty` for terminal bar charts, `--json` for machine-readable
output (an array of result objects), or `-o results.json` to write to a
file. Full reference: [docs/cli.md](docs/cli.md).

## Supported input formats

| Format      | Value (`-v`)  | Notes                                        |
| ----------- | ------------- | -------------------------------------------- |
| 23andMe     | `23andme`     | whitespace, 4 columns                        |
| AncestryDNA | `ancestry`    | whitespace, 5 columns                        |
| FTDNA       | `ftdna`       | CSV, 4 columns                               |
| FTDNA alt.  | `ftdna2`      | CSV, 5 columns                               |
| WeGene      | `wegene`      | same layout as 23andMe                       |
| MyHeritage  | `myheritage`  | same layout as FTDNA                         |
| LivingDNA   | `livingdna`   | same layout as 23andMe                       |
| VCF / gVCF  | `vcf`/`gvcf`  | **auto-detected** — rsID or position matching, hg19+hg38 |

## Available models (`-m`)

27 models ship with data embedded — the Dodecad family (`K7b`, `K12b`,
`globe13`, `globe10`, `world9`, `Eurasia7`, `Africa9`, `weac2`), Eurogenes
(`K36`, `EUtest13`, `Jtest14`, `K7AMI`, `K8AMI`), `HarappaWorld`,
fine-grained global calculators (`E11`, `MDLPK27`, `K47`, `MichalK25`), the
Tolan family (`K7M1`, `K13M2`, `K14M1`, `K18M4`, `K25R1`), regional ones
(`TurkicK11`, `KurdishK10`, `puntDNAL`, `AncientNearEast13`) — 35 to 47
populations, 27k to 188k SNPs each.

```bash
./admix -f genome.txt -m K12b,K36        # several models in one run
./admix -f genome.txt                     # everything (slow, verbose)
```

Full catalog with SNP counts, population components and links to the
original calculators: [docs/models.md](docs/models.md).

## API server

```bash
go build -o server ./cmd/server
./server                # listens on 127.0.0.1:2364
./server -port 9000     # custom port
```

```bash
curl -F "raw_data=@genome.txt" -F "model=K12b" -F "vendor=23andme" \
     http://127.0.0.1:2364/calculate
```

`GET /models` lists every registered model, `GET /health` is a liveness
probe. Full reference with request and response examples:
[docs/api.md](docs/api.md).

## Project layout

```
cmd/admix/         CLI entry point
cmd/server/        HTTP API entry point
cmd/snpindex/      coordinate-index generator (dev tool)
pkg/parser/        vendor text parsers + VCF/gVCF parser
pkg/calculator/    likelihood model + gradient descent optimizer
pkg/models/        model registry, loaders, coordinate index, embedded data
pkg/result/        shared result shape (CLI, API, renderer)
pkg/render/        terminal bar-chart renderer (-pretty)
build.sh           quick local cross-build (releases use GoReleaser)
```

## Development

```bash
make test       # unit + integration tests
make race       # tests with the race detector
make bench      # parser benchmarks
make lint       # gofmt + golangci-lint
```

CI runs the test matrix (Linux/Windows/macOS, race detector), vet and
golangci-lint on every push/PR; pushing a `v*` tag publishes a draft
release with binaries for linux/windows/macos (amd64/arm64) via GoReleaser.

## Documentation

- [docs/architecture.md](docs/architecture.md) — pipeline, packages,
  optimization algorithm, coordinate index
- [docs/cli.md](docs/cli.md) — CLI reference, all input formats,
  VCF/gVCF behavior and coverage, troubleshooting
- [docs/api.md](docs/api.md) — HTTP API reference
- [docs/models.md](docs/models.md) — model catalog, data file formats,
  adding new models

## Credits & data sources

This project started as a Go port of
[stevenliuyi/admix](https://github.com/stevenliuyi/admix) (Python,
GPL-3.0): the model catalog, the population registry, the embedded data
file formats and the likelihood formulation follow that project, and the
model data files are taken from its compilation. Everything else — the Go
implementation, VCF/gVCF parsing with positional matching on hg19/hg38,
the coordinate index, the HTTP API and the terminal renderer — is original
to this repository.

The admixture calculators themselves are community creations published on
blogs and forums (Dodecad, Eurogenes/BGA101, HarappaDNA, MDLP, Anthrogenica
and others). **All calculator files remain the property of their original
authors**; per-model source links are listed in
[docs/models.md](docs/models.md#model-catalog). If you are one of the
original authors and want your calculator removed from this distribution,
please open an issue.

This program is GPL-3.0 like the project it derives from; see
[LICENSE](LICENSE).

## Limitations

- Estimates are model-dependent interpretations, not facts; different
  calculators decompose the same genome very differently. Treat results as
  entertainment and a starting point for genealogical research — not as
  medical, legal or identity evidence.
- No phasing, no strand alignment: your genotype is matched as-is against
  each model's allele definitions.
- Exome and other partial VCFs cover only the sites they assayed. A
  whole-genome variant VCF is filled in as homozygous-reference at the
  model SNPs it does not list; pass `-homref off` to keep those sites
  missing. The CLI reports the matched/total ratio either way. See
  [docs/cli.md](docs/cli.md#vcfgvcf-input).

## License

[GPL-3.0](LICENSE)
