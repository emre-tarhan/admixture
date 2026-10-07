# HTTP API reference

The `server` binary exposes the same calculation pipeline as the CLI over
HTTP. It is designed for local use — a GUI frontend on the same machine, a
small intranet tool, etc.

## Running

```bash
go build -o server ./cmd/server
./server                     # 127.0.0.1:2364
./server -port 8080          # custom port
PORT=8080 ./server           # via environment (override with -port)
./server -addr 0.0.0.0      # expose — only behind an authenticating proxy
```

- Default port **2364**, address **127.0.0.1** (not reachable from other
  hosts unless `-addr 0.0.0.0` is passed).
- No authentication, no TLS, no rate limiting — do not put it directly on a
  public network.

## Endpoints

### `GET /health`

Liveness probe: `{"status":"ok"}`.

### `GET /models`

Lists every registered model name (including ones without embedded data —
see the note in [models.md](models.md#registered-but-without-data)).

```bash
curl http://127.0.0.1:2364/models
```

```json
{"models":["K7b","K12b","globe13","globe10","world9","Eurasia7","Africa9","weac2","E11","K36","EUtest13","Jtest14","HarappaWorld","TurkicK11","KurdishK10","AncientNearEast13","K7AMI","K8AMI","MDLPK27","puntDNAL","K47","K7M1","K13M2","K14M1","K18M4","K25R1","MichalK25"]}
```

### `POST /calculate`

Runs the admixture calculation on an uploaded raw data file.

Request: `multipart/form-data` with:

| Field      | Required | Default    | Description                                             |
| ---------- | -------- | ---------- | ------------------------------------------------------- |
| `raw_data` | yes      | —          | The raw genome file (everything the CLI accepts: vendor text formats, VCF, gVCF). |
| `model`    | no       | `K7b`      | Model name, or comma-separated names (e.g. `K12b,K36`). |
| `vendor`   | no       | `23andme`  | Vendor format key — same values as the CLI `-v` flag ([cli.md](cli.md#vendor-file-formats)). VCF/gVCF uploads are auto-detected regardless of this field. |
| `build`    | no       | `auto`     | Reference build for VCF position matching (`auto`, `hg19`, `hg38`). |

Example:

```bash
curl -F "raw_data=@genome.vcf.gz" \
     -F "model=K12b,K36" \
     -F "build=auto" \
     http://127.0.0.1:2364/calculate
```

Response: `application/json` — an **array** of result objects, one per
model. Proportions are raw shares (0–1); `matched_snps`/`total_snps` let
the client judge input coverage (see [cli.md](cli.md#vcfgvcf-input)).

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

#### Error semantics

- Non-POST to `/calculate` → `405 Method not allowed`.
- Malformed multipart body → `400 "Unable to parse form"`.
- Missing `raw_data` part → `400 "Missing raw_data file"`.
- Request naming only unknown/undataed models →
  `400 "No requested model is available"`.
- Unparseable genome file (wrong vendor format etc.) →
  `400 "Error parsing data: ..."`.
- **Unknown or undataed model names are silently skipped** — a request for
  `model=K12b,Globe13` (Globe13 has no data embedded) returns `200` with
  only K12b in the array. Clients should check that every requested model
  is present in the response.
- The optimizer tolerance is fixed at `1e-4` for API requests (the CLI's
  `-t` flag has no API equivalent).

## Limits and lifecycle

- Uploads are capped at **10 MB** (`ParseMultipartForm(10 << 20)`); larger
  bodies fail with `400`. Note: whole-genome VCFs are usually far bigger
  than 10 MB — the CLI is the right tool for those. Vendor raw text files
  (~100 MB uncompressed, ~15–30 MB gzipped) typically fit.
- Each upload is written to a temporary file (the parsers are file-based)
  and **deleted when the request handler returns**. Genome data is not
  persisted, logged or sent anywhere by the server itself.
- Calculation is synchronous — one request holds one goroutine for the full
  optimizer run. Concurrent requests run concurrently; there is no queue or
  timeout of its own (set one on your HTTP client / reverse proxy).
