# Models

An admixture **model** defines a set of reference populations and, for every
SNP it covers, the minor/major allele pair and the frequency of the major
allele in each population. The optimizer (see
[architecture.md](architecture.md#the-estimation-problem)) finds the
mixture of those populations that best explains your genotypes.

Most bundled models are well-known community calculators in the tradition
of the Dodecad, Eurogenes, Harappa and MDLP projects — the same models many
people know from GEDmatch-style analysis.

## Model catalog

All 27 registered models ship with data embedded in the binary. Frequency
matrices are stored gzip-compressed to keep the binary at ~205 MB
(~190 MB in release archives):

| Model               | SNPs       | Pops | Calculator                      | Source |
| ------------------- | ---------- | ---- | ------------------------------- | ------ |
| `K7b`               | 166,770    | 7    | Dodecad K7b                     | [Dodecad](http://dodecad.blogspot.com/2012/01/k12b-and-k7b-calculators.html) |
| `K12b`              | 166,770    | 12   | Dodecad K12b                    | [Dodecad](http://dodecad.blogspot.com/2012/01/k12b-and-k7b-calculators.html) |
| `globe13`           | 166,255    | 13   | Dodecad globe13                 | [Dodecad](http://dodecad.blogspot.com/2012/10/globe13-calculator.html) |
| `globe10`           | 166,255    | 10   | Dodecad globe10                 | [Dodecad](http://dodecad.blogspot.com/2012/10/globe10-calculator.html) |
| `world9`            | 170,822    | 9    | Dodecad world9                  | [Dodecad](http://dodecad.blogspot.com/2011/12/world9-calculator.html) |
| `Eurasia7`          | 164,990    | 7    | Dodecad Eurasia7                | [Dodecad](http://dodecad.blogspot.com/2011/10/eurasia7-calculator.html) |
| `Africa9`           | 27,187     | 9    | Dodecad Africa9                 | [Dodecad](http://dodecad.blogspot.com/2011/09/africa9-calculator.html) |
| `weac2`             | 169,131    | 7    | Dodecad weac2 (West Eurasian cline) | [Dodecad](http://dodecad.blogspot.com/2012/06/weac2-calculator.html) |
| `E11`               | 170,288    | 11   | E11                             | [WeGene/Ranhaer](http://www.ranhaer.com/thread-32241-1-1.html) |
| `K36`               | 165,688    | 36   | Eurogenes K36                   | [BGA101](http://bga101.blogspot.com/2013/03/eurogenes-k36-at-gedmatch.html) |
| `EUtest13`          | 159,620    | 13   | Eurogenes EUtest K13            | [BGA101](https://bga101.blogspot.com/2013/11/updated-eurogenes-k13-at-gedmatch.html) |
| `Jtest14`           | 159,620    | 14   | Eurogenes Jtest K14 (Ashkenazim)| [BGA101](http://bga101.blogspot.com/2012/09/eurogenes-ashkenazim-ancestry-test-files.html) |
| `HarappaWorld`      | 188,173    | 16   | HarappaWorld                    | [HarappaDNA](http://www.harappadna.org/2012/05/harappaworld-admixture/) |
| `TurkicK11`         | 112,117    | 11   | Turkic K11                      | [Anthrogenica](http://www.anthrogenica.com/showthread.php?8817-Turkic-K11-Admixture-Calculator) |
| `KurdishK10`        | 92,632     | 10   | Kurdish K10                     | [Anthrogenica](http://www.anthrogenica.com/showthread.php?8571-K10-Kurdish-Calculator-Version-1/page6) |
| `AncientNearEast13` | 35,648     | 13   | Ancient Near East K13           | [Anthrogenica](http://www.anthrogenica.com/showthread.php?8193-ancient-DNA-in-the-Gedrosia-Near-East-Neolithic-K13) |
| `K7AMI`             | 101,647    | 7    | Eurogenes K7 AMI                | [Anthrogenica](http://www.anthrogenica.com/showthread.php?4548-Upcoming-DIY-Eurogenes-K7-amp-K8-Calculator-amp-Oracles-for-tracking-E-Asian-amp-ASI) |
| `K8AMI`             | 101,647    | 8    | Eurogenes K8 AMI                | [Anthrogenica](http://www.anthrogenica.com/showthread.php?4548-Upcoming-DIY-Eurogenes-K7-amp-K8-Calculator-amp-Oracles-for-tracking-E-Asian-amp-ASI) |
| `MDLPK27`           | 118,536    | 27   | MDLP K27                        | [Anthrogenica](http://www.anthrogenica.com/showthread.php?4557-Post-MDLP-K27-Results) |
| `puntDNAL`          | 38,496     | 12   | puntDNAL K12 Ancient World      | [Anthrogenica](http://www.anthrogenica.com/showthread.php?8034-PuntDNAL-K12-Ancient-World-Results) |
| `K47`               | 76,267     | 47   | LM Genetics K47                 | [Anthrogenica](https://anthrogenica.com/showthread.php?12788-New-K30-K47-World-Calculator) |
| `K7M1`              | 165,688    | 7    | Tolan K7M1                      | [Gen3553](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K13M2`             | 165,688    | 13   | Tolan K13M2                     | [Gen3553](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K14M1`             | 165,688    | 14   | Tolan K14M1                     | [Gen3553](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K18M4`             | 165,688    | 18   | Tolan K18M4                     | [Gen3553](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `K25R1`             | 165,688    | 25   | Tolan K25R1                     | [Gen3553](http://gen3553.pagesperso-orange.fr/ADN/Calc.htm) |
| `MichalK25`         | 76,261     | 25   | Michal World K25                | [Anthrogenica](https://anthrogenica.com/showthread.php?13359-Michal-s-World-K25-calculator) |

Population components for every model (in registry order — the column
order of the `.F` files) are defined in
[`pkg/models/registry.go`](../pkg/models/registry.go).

Note: Eurogenes K13 data files exist in the upstream compilation, but the
model is deliberately not exposed here or there — its population order
could not be verified, and a wrong column mapping silently produces wrong
results.

These are community calculators originally published on blogs and forums
(Dodecad, Eurogenes/BGA101, HarappaDNA, MDLP, Anthrogenica threads...).
All calculator files remain the property of their original authors; see
[Credits](../README.md#credits--data-sources).

## Data file formats

Model data lives in `pkg/models/data/` and is embedded into the binary at
compile time. Two text files per model, whitespace-separated:

**`<Model>.alleles`** — one line per SNP:

```
rsID minor major
rs1000000 T C
rs10000023 G T
```

**`<Model>.<N>.F.gz`** — `<N>` is the number of populations (e.g.
`K12b.12.F.gz`). Same content as a plain `.F` file, stored **gzip-compressed**
and decompressed transparently on load. One row per SNP, **in the same order
as the `.alleles` file**, with `N` whitespace-separated frequencies, one per
population:

```
0.705458 0.861519 0.933294 ...   ← row for rs1000000
0.544086 0.408918 0.529086 ...   ← row for rs10000023
```

The frequency is that of the **major** allele — the allele in the third
column of the `.alleles` file. (Note: the matrix values routinely exceed
0.5, which rules out a minor-allele interpretation.) Column *k* corresponds
to the population at index *k* of that model's list in `GetPopulations`.

`LoadModel` validates on load: every frequency row must have exactly `N`
columns, and the SNP count must equal the row count — otherwise the model
is rejected with an error.

## Coordinate index (`positions.tsv`)

Vendor text formats are matched by rsID, which is build-independent. VCF
and gVCF files, however, often carry no rsIDs — matching then happens by
chromosomal position, which requires knowing where each model SNP sits on
each reference build.

`data/positions.tsv` provides that mapping, one line per SNP:

```
rsID  chrom  pos19  pos38  ref
rs10000023  4  95733906  94812755  G
```

- `pos19` / `pos38` — 1-based position on hg19 / hg38 (`.` when unknown)
- `ref` — expected reference allele on the primary assembly (forward
  strand), used to resolve homozygous-reference calls inside gVCF
  reference blocks and, for whole-genome variant VCFs, at model SNPs the
  file does not list. `.` means the allele is unknown and no
  homozygous-reference call is invented for that SNP.

The file carries a primary-chromosome coordinate on both builds for
**332,657 of the 333,056** unique model SNPs. It is generated by
[`cmd/snpindex`](../cmd/snpindex/main.go)
from the UCSC `snp151Common` tables (dbSNP build 151) for both assemblies.
When an rsID has several rows, a primary chromosome (1–22, X, Y, M) is
kept and a later alt-contig row does not replace it. Rows that had
already been stored on an alt contig were moved back onto the primary
assembly by lifting their hg19 position with UCSC `hg19ToHg38`; the
reference allele on those rows is `.`, because the alt-contig base is
often not the primary reference.

```bash
go run ./cmd/snpindex \
  -alleles 'pkg/models/data/*.alleles' \
  -hg19 snp151Common_hg19.txt.gz \
  -hg38 snp151Common_hg38.txt.gz \
  -out pkg/models/data/positions.tsv
```

## Adding a model

1. Append the population list to the `switch` in
   [`pkg/models/registry.go`](../pkg/models/registry.go)
   (`GetPopulations`). Keep the order fixed — it defines the `.F` column
   order.
2. Put `<Model>.alleles` into `pkg/models/data/` and gzip the frequency
   matrix to `<Model>.<N>.F.gz` (`N` = number of populations you
   registered) — the loader opens only the gzip-compressed name.
3. Re-run `cmd/snpindex` (command above) so the new SNPs get coordinates.
4. Rebuild — `//go:embed data/*` picks up new files automatically.

Practical notes:

- Source data for most classic calculators circulates in the community as
  exactly these two file types (many GEDmatch-derived tools use the same
  format); if you have a `<Model>.txt`/`<Model>.<N>.F` pair from another
  tool, verify the column semantics before reusing — the frequency matrix
  must be the major-allele frequency in registry column order.
- Models with more SNPs are slower (cost per iteration is
  `O(#SNPs × #populations)`) but generally more stable; models below a few
  thousand SNPs get noisy.
- After adding, verify with:
  `go run ./cmd/admix -f <rawfile> -m <Model> -v 23andme`.
