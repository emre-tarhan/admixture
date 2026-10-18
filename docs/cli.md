# CLI reference

`admix` reads a raw autosomal data file and prints admixture proportions for
one or more models.

## Synopsis

```
admix -f <file> [-m <models>] [-v <vendor>] [-b <build>] [-homref <auto|on|off>] [-o <file>] [--json | -pretty] [-t <tolerance>]
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
| `-homref`| `auto`    | How to treat model SNPs a variant-only VCF does not list. `auto` calls them homozygous-reference on a whole-genome VCF (≥1,000,000 records on ≥20 autosomes). `on` does that for any variant-only VCF. `off` leaves them missing. gVCFs are never filled in, because their reference blocks already say where the sample is homozygous-reference. |
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

The matched/total ratio is the single most important quality signal: chips,
vendor exports, gVCFs and whole-genome variant VCFs (with the default
`-homref auto`) match ~95%+. Exomes and other partial VCFs stay near the
share of the model that the file actually assayed. Low ratios mean the
optimizer is working with less information — see
[the VCF notes below](#vcfgvcf-input).

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
  span. No-call spans (`./.`) are skipped. A file that contains any
  reference block is not otherwise filled in: a gap between blocks stays
  "not genotyped".
- **Whole-genome variant VCFs** — a variant-only callset has no record
  where the sample matches the reference. `-homref auto` (the default)
  treats that as homozygous-reference when the file looks like a whole
  genome: at least 1,000,000 data records and variants on at least 20
  autosomes. The genotype is the reference allele from the coordinate
  index. This applies only on primary chromosomes where at least one
  model SNP was actually seen, so a wrong `-b` value cannot turn the
  whole genome into reference calls. `-homref on` applies the same rule
  to smaller files (an exome, a chromosome slice). `-homref off` keeps
  every absent site missing.
- Filtered records (`FILTER` other than `PASS`/`.`), indels, MNPs,
  symbolic alleles and partial no-calls are skipped, and they are not
  replaced by a homozygous-reference guess.
- Multi-sample VCFs: the **first** sample column is used.

Measured on one DRAGEN hg38 whole-genome VCF (5.1M records, no rsIDs)
against a 23andMe v3 file built from that same genome. K12b:

| Input                                      | Matched SNPs       |
| ------------------------------------------ | ------------------ |
| WGS variant VCF, `-homref auto` (default)  | 164,974 / 166,770  |
| Same file, `-homref off`                   | 86,998 / 166,770   |
| 23andMe v3 from the same genome            | 166,769 / 166,770  |
| 30× gVCF, hg38 (reference blocks)          | 163,648 / 166,770  |
| 1000 Genomes chr22 slice (rsIDs, one chrom)| 2,546 (chr22 only) |

On this pair the default VCF run and the 23andMe run describe the same
person: every K12b component agrees within 0.2 percentage points, and
every K36 component within 0.4. `-homref off` on the same VCF moves
K12b Caucasus from 43.0% to 35.2% and invents about 10% African that
the chip file does not have. A single-chromosome slice stays under the
size rule, so `auto` does not fill the other chromosomes. The gVCF row
is an earlier 30× sample whose reference blocks already supplied the
homozygous-reference calls; that path is unchanged.

**Coverage.** With `-homref auto`, a whole-genome variant VCF and a gVCF
both resolve homozygous-reference model SNPs, so coverage lands in the
same range as a vendor export. The assumption is only as good as the
callset: a site filtered out, or sitting in a region the caller never
looked at, must not be filled in. Exomes and chromosome slices therefore
stay below the size rule and keep absent sites missing — pass `-homref on`
only when you know the file is a complete variant list for the regions
you care about. `-homref off` is the pre-0.1.1 behavior, and it biases
the fit, because the dropped sites are mostly homozygous for the common
allele.

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
