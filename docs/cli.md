# CLI reference

`admix` reads a raw autosomal data file and prints admixture proportions for
one or more models.

## Synopsis

```
admix -f <file> [-m <models>] [-v <vendor>] [-b <build>] [-o <file>] [--json | -pretty] [-t <tolerance>]
admix -models
admix --help
```

## Flags

| Flag     | Default   | Description                                                        |
| -------- | --------- | ------------------------------------------------------------------ |
| `-f`     | *(none)*  | **Required.** Path to the raw genome data file. `.vcf`, `.vcf.gz`, `.vcf.bgz` and gzip'd text files are all fine. |
| `-models`| off       | List all available models (with population counts) and exit.       |
| `-m`     | all       | Comma-separated model names, e.g. `K12b,K36`. Omit to run every registered model. |
| `-v`     | `23andme` | Vendor format of the input file. One of `23andme`, `ancestry`, `ftdna`, `ftdna2`, `wegene`, `myheritage`, `livingdna`, `vcf`, `gvcf`. VCF/gVCF input is auto-detected by content, so the flag only matters for text files. |
| `-b`     | `auto`    | Reference build for VCF position matching: `auto`, `hg19` or `hg38`. `auto` detects the build from `##contig` lengths and only matters when the file has no rsIDs. |
| `-o`     | *(stdout)*| Write results to this file instead of stdout.                      |
| `--json` | off       | Emit JSON (an array of result objects) instead of formatted text.  |
| `-pretty`| off       | Render results as terminal bar charts (see below). Mutually exclusive with `--json`. |
| `-t`     | `1e-4`    | Optimizer convergence tolerance (loss improvement threshold). Higher = faster, less precise. |

## Model selection

All 27 registered models ship with data (catalog with sources:
[models.md](models.md)). Highlights:

- `K12b`, `K7b`, `globe13`, `globe10`, `world9`, `Eurasia7`, `Africa9`,
  `weac2` — Dodecad calculators
- `K36`, `EUtest13`, `Jtest14` — Eurogenes; `HarappaWorld` — HarappaDNA
- `E11` — East Asia focus; `MDLPK27`, `K47`, `MichalK25` — fine-grained
  global; `TurkicK11`, `KurdishK10`, `AncientNearEast13` — regional
- `K7AMI`, `K8AMI`, `K7M1`, `K13M2`, `K14M1`, `K18M4`, `K25R1`,
  `puntDNAL` — ancient/ancient-adjacent components

- `-m K12b` — single model.
- `-m K12b,K36,AncientNearEast13` — several models in one run; results are
  printed grouped per model.
- No `-m` — every registered model runs (this takes a while and produces a
  lot of output; prefer naming the models you want).

## Output formats

### Text (default)

Populations sorted by percentage, descending; contributions below 0.01% are
omitted; a blank line separates models. Every model block ends with a SNP
match statistic — how many of the model's SNPs were found in your input:

```
Parsing raw data from genome.txt (23andme)...
Loaded 610526 genotypes
Calculating K12b...
Model: K12b
Caucasus: 30.92%
North European: 26.61%
Atlantic Med: 19.44%
(121,334 of 166,770 SNPs matched)
```

The matched/total ratio is the single most important quality signal: chips
and vendor exports match ~95%+, plain WGS VCFs (variant sites only) land
around 50%, gVCFs with reference blocks ~98%. Low ratios mean the
optimizer is working with less information — see
[vcf.md notes below](#vcfgvcf-input).

### Bar charts (`-pretty`)

```
 K12b · 163,648 / 166,770 SNPs matched
────────────────────────────────────────────────────────────────
  Caucasus        █████████████████████████░  25.25%
  South Asian     ██████████████████████░░░░  22.24%
  Gedrosia        █████████████████░░░░░░░░░  16.74%
  ...
```

Bars are relative to the largest population in that model. Colors are
enabled automatically on terminals and disabled on pipes/redirects and
when `NO_COLOR` is set.

### JSON (`--json`)

An array of result objects — proportions are raw shares (0–1), zero values
included:

```json
[
  {
    "model": "K12b",
    "proportions": {
      "Caucasus": 0.2525,
      "South Asian": 0.2224,
      "Gedrosia": 0.1674
    },
    "matched_snps": 163648,
    "total_snps": 166770
  }
]
```

With `-o results.json` the output is written to the file and the terminal
shows only progress lines.

## Vendor file formats

Which `-v` value matches which download, and what the parser expects:

| `-v` value   | Service                    | Expected layout                                            |
| ------------ | -------------------------- | ---------------------------------------------------------- |
| `23andme`    | 23andMe raw data           | `rsID chromosome position genotype`, whitespace-separated   |
| `wegene`     | WeGene                     | same layout as 23andMe                                     |
| `livingdna`  | LivingDNA                  | same layout as 23andMe                                     |
| `ancestry`   | AncestryDNA                | `rsID chromosome position allele1 allele2`, whitespace-separated; header line skipped |
| `ftdna`      | FamilyTreeDNA / MyHeritage | `RSID,CHROM,POS,RESULT` CSV; quoted fields tolerated       |
| `myheritage` | MyHeritage                 | same layout as `ftdna`                                     |
| `ftdna2`     | FamilyTreeDNA (alternate export) | `rsID,chrom,pos,allele1,allele2` CSV, 5 columns      |
| `vcf`/`gvcf` | any VCF/gVCF               | auto-detected by content — see next section                |

Notes:

- Comment lines starting with `#` and blank lines are ignored for all
  formats.
- Separators are auto-fallback: the FTDNA parsers try CSV first and fall
  back to whitespace-splitting when a line has fewer than the expected
  number of comma-separated fields.
- The genotype is taken as-is (first and last character of the genotype
  string); there is no strand flipping or position alignment against the
  model files for text formats.

## VCF/gVCF input

VCF and gVCF files (plain, `.gz`/bgzip-compressed) are auto-detected and
parsed with position-aware matching:

- **rsID matching** — records with a dbSNP ID (`ID` column) are matched by
  rsID, regardless of build. Most published VCFs (1000 Genomes, dbSNP-based
  callsets, vendor VCF exports) work this way.
- **Position matching** — WGS callsets often have an empty `ID` column.
  For those, records are matched positionally against each model SNP's
  coordinates on **both hg19 and hg38** (embedded in the binary; generated
  from UCSC dbSNP b151 via `cmd/snpindex`). The reference build is
  auto-detected from the header's `##contig` lengths; pass `-b hg19` /
  `-b hg38` to override (only needed for header-less files).
- **gVCF reference blocks** — hom-ref spans (`ALT=<*>`/`<NON_REF>`,
  `END=`) resolve homozygous-reference calls for model SNPs inside the
  span, roughly doubling matched SNPs versus a variant-only VCF. No-call
  spans (`./.`) are skipped.
- Filtered records (`FILTER` other than `PASS`/`.`), indels, MNPs,
  symbolic alleles and partial no-calls are skipped.
- Multi-sample VCFs: the **first** sample column is used.

Measured behavior (K12b, 8-core VPS):

| Input                                            | Matched SNPs      |
| ------------------------------------------------ | ----------------- |
| WGS variant VCF, hg38, no rsIDs (~4.0M records)  | 85,257 / 166,770  |
| Same sample lifted to hg19                       | 85,256            |
| 30× gVCF, hg38 (with reference blocks)           | 163,648 / 166,770 |
| Same sample, hg19 variant-only export            | 85,227            |
| 1000 Genomes chr22 slice, 2504 samples (rsIDs)   | 2,546 (chr22 only)|

The hg38 and hg19 files of the same sample produce near-identical
proportions (largest observed difference 0.06 percentage points), which
confirms position matching works on both builds and that `bcftools
liftOver` output is handled cleanly.

**Coverage caveat.** A variant-only VCF contains no record for positions
where the sample is homozygous-reference, so roughly half of a model's
SNPs are missing — and the missingness is population-correlated (common
alleles in your ancestry are the ones you're homozygous-reference for).
Treat VCF-based results as approximations; a gVCF or a vendor raw export
is always the better input. gVCF input is the most complete (hom-ref
resolvable), followed by vendor text files (full chip coverage), then
variant-only VCFs.

## Exit codes

- `0` — success (even if some selected models failed and reported errors).
- `1` — missing `-f`, unreadable input file, a parse error, or no
  requested model could be loaded.

## Troubleshooting

**"no raw data file given (-f)"** — the flag is required.

**"Loaded 0 genotypes"** — the file parsed but nothing matched the
expected column layout: wrong `-v` for the file, a header-only or empty
file, or a format the parser doesn't recognize (e.g. BAM). Open the file
and compare against the layouts above.

**"VCF, build undetermined — rsID matching only"** — the VCF header has no
`##contig` lines with recognizable lengths. Pass `-b hg19` or `-b hg38` to
enable position matching (rsID matching still worked).

**Very low matched-SNP ratio on a VCF** — likely one of: (a) the file is a
variant-only VCF (expected — see coverage caveat above), (b) wrong build
with no rsIDs (force with `-b`), or (c) the file is not hg19/hg38 at all
(e.g. a different species or a non-standard assembly).

**"Error loading model X: model X not found"** — `X` is not a registered
model name. See the catalog ([models.md](models.md#model-catalog)) for
valid names.

**Results look implausible (e.g. one population near 100%)** — almost
always a vendor-format mismatch: a wrong `-v` value makes most genotypes
parse as missing or mismatched, and the optimizer then fits noise. Check
the matched-SNP ratio first, then `-v`.
