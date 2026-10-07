# Architecture

This document describes how admixture turns a raw genotype file into
population proportions, and how the code is organized.

## Pipeline

```
raw data file (vendor text format, or VCF/gVCF)
        │
        ▼
┌─────────────────┐   lookup key → genotype
│ pkg/parser      │   keys: "rs3094315" (rsID) and/or
│                 │   "hg38:1:752566" (build:chrom:pos)
└─────────────────┘
        │
        ▼
┌─────────────────┐   per-SNP major/minor allele counts
│ pkg/calculator  │   gMajor[i], gMinor[i] ∈ {0, 1, 2}
│  genotypeMatch  │
└─────────────────┘
        │
        ▼
┌─────────────────┐   mixture vector q (|q| = #populations,
│ pkg/calculator  │   Σq = 1, q ≥ 0) + match statistics
│  optimizer      │
└─────────────────┘
        │
        ▼
   text / bar chart / JSON output
```

## Packages

### `pkg/parser`

Two parsing paths behind one entry point:

- `ParseFile(filename, vendor, wanted, build)` sniffs the file content:
  files whose first line is `##fileformat=VCF` (optionally gzipped/bgziped)
  go to the VCF parser, everything else to the vendor text parsers.
- Text parsers (one per vendor layout) return `map[rsID]genotype`. See
  [cli.md](cli.md#vendor-file-formats) for the exact layouts.
- `ParseVCF` handles both plain VCF and gVCF. Matching is done against a
  `models.Wanted` key set derived from the selected models:
  - Records carrying a dbSNP ID emit their **rsID key** (build-agnostic).
  - When the reference build is known, records also emit their
    **coordinate key** (`hgNN:chrom:pos`) — this covers WGS callsets that
    carry no rsID annotations at all.
  - The build is auto-detected from `##contig` length declarations
    (hg19 and hg38 chromosome lengths never overlap) and can be forced
    with `-b hg19|hg38`.
  - gVCF **reference blocks** (ALT `<*>`/`<NON_REF>`, INFO `END=`) assert
    the sample is homozygous-reference across their span: every wanted
    coordinate inside the span whose expected reference allele is known is
    emitted as a hom-ref genotype. No-call blocks (`./.`) are skipped.
  - Indels, MNPs, symbolic alleles, non-PASS records and partial no-calls
    (`0/.`, `./.`) are skipped.

### `pkg/models`

Holds model definitions, the coordinate index and the embedded data:

- `registry.go` — hard-coded list of populations per model
  (`GetPopulations`) and the full model list (`ListModels`).
- `embed.go` — `//go:embed data/*` pulls every file in `data/` into the
  binary at compile time.
- `coords.go` — the coordinate index (`data/positions.tsv`): rsID →
  chromosomal position on hg19 and hg38 plus the expected reference allele.
  `BuildWanted` unions the lookup keys (rsIDs and coordinate keys) of the
  selected models into the set the VCF parser extracts; range queries over
  that set power gVCF reference-block expansion.
- `loader.go` — `LoadModel(name)` reads two files from the embedded FS and
  validates them against each other:
  - `data/<Model>.alleles` — one SNP per line: `rsID minor major`
  - `data/<Model>.<N>.F` — one row per SNP (same order) with `N`
    whitespace-separated frequencies, one per population. The frequency is
    that of the **major** allele (the third column of `.alleles`) —
    verified empirically: values routinely exceed 0.5, which is impossible
    for minor-allele frequencies.

  File formats are specified in [models.md](models.md).

### `pkg/calculator`

`calculator.Calculate(model, rawData, tolerance)` runs the estimation. It is
pure math — no I/O — and returns the mixture plus `Stats` (how many model
SNPs were actually found in the input, the key quality signal for VCF
input).

### `pkg/result` / `pkg/render`

`result.Result` is the shared output shape (model, proportions, SNP match
counts) used by CLI, API and renderer. `render.Pretty` draws the terminal
bar charts (`-pretty`) with plain ANSI codes — no color dependency, and
color is disabled automatically on non-TTY output and when `NO_COLOR` is
set.

## The estimation problem

For each model SNP *i* with per-population major-allele frequencies
`F[i][k]` (population *k*), the count of major alleles observed in the
individual is modeled as Binomial(2, `p[i]`) with

```
p[i] = Σₖ F[i][k] · q[k]
```

Given observed counts `gMajor[i]` (major-allele copies, 0–2) and
`gMinor[i]` (minor-allele copies, 0–2), the optimizer minimizes the negative
log-likelihood

```
L(q) = −Σᵢ [ gMajor[i]·log(p[i]) + gMinor[i]·log(1−p[i]) ]
```

over the probability simplex (`q ≥ 0`, `Σq = 1`). The gradient is

```
∂L/∂q[k] = −Σᵢ [ gMajor[i]/p[i] − gMinor[i]/(1−p[i]) ] · F[i][k]
```

`p[i]` is clamped to `[1e-9, 1−1e-9]` for numerical stability, and SNPs
with no observed data (`gMajor = gMinor = 0`, i.e. the rsID/coordinate was
not found in the input) are skipped entirely.

## The optimizer

**Projected gradient descent** with an adaptive step size:

1. Initialize `q` uniformly (`1/|populations|` each).
2. Take a gradient step, then project back onto the simplex using the
   standard sorting algorithm (Wang & Carreira-Perpiñán 2014,
   <https://eng.ucmerced.edu/people/wwang5/papers/SimplexProj.pdf>).
3. If the loss improved: accept, multiply the learning rate by **1.1**.
   Otherwise: reject (keep old `q`), multiply by **0.5**.
4. Stop when the loss improvement drops below `-t` (default `1e-4`), the
   learning rate falls below `1e-10`, or after **5000** iterations.

Because the objective is convex and the feasible set is the simplex, this
converges reliably from the uniform start. Cost per iteration is
`O(#SNPs × #populations)`; a full K12b fit takes well under a second. The
per-SNP loop is currently sequential — parallelizing across SNPs is the
obvious next optimization if larger models are added.

## Binaries

| Binary     | Entry point      | Role                                                          |
| ---------- | ---------------- | ------------------------------------------------------------- |
| `admix`    | `cmd/admix`      | CLI: flags, prints/writes text, bar charts or JSON results     |
| `server`   | `cmd/server`     | HTTP API: multipart upload → temp file → same pipeline         |
| `snpindex` | `cmd/snpindex`   | Dev tool: regenerates `data/positions.tsv` from UCSC dbSNP tables |

The server writes each upload to an `os.CreateTemp` file (the parsers are
file-based) and deletes it when the request finishes. See [api.md](api.md).

## Input reference builds (hg19 vs hg38)

Vendor text formats are matched by rsID only — the reference build does not
matter for them. VCF/gVCF input uses positions when rsIDs are missing, so
the build matters there; both builds are supported simultaneously because
every model SNP carries both coordinates in the embedded index, and the
build used is auto-detected per file (or forced via `-b`). A lifted-over
VCF (e.g. `bcftools liftOver` hg38 → hg19) therefore produces essentially
identical results; see [cli.md](cli.md#vcfgvcf-input) for measured
comparisons.

## Design notes

- **Embedding** — `go:embed` keeps the whole tool a single static file
  (models + coordinate index inside the binary).
- **No dependencies** — `go.mod` has zero requires; builds offline.
- **Error strategy** — the CLI reports per-model load errors and continues
  with the remaining models; the API server skips failing models (see
  [api.md](api.md#error-semantics)).
- **Registered vs. available models** — every one of the 27 registered
  models ships with data; loading a name that is not in `ListModels()`
  returns "model not found".
