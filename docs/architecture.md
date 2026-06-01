# Architecture

This document describes how the platform's pieces fit together and the
contracts between them.

## Components

| Component | Tech | Port | Role |
|-----------|------|------|------|
| `frontend` | Next.js 14, TS, Tailwind | 3000 | Product UI (auth, dashboard, analyses, results, reports, admin) |
| `gateway` | FastAPI, SQLAlchemy, JWT | 8080 | Single entry point: auth, persistence, orchestration, fusion, reports |
| `genomics-validation` | FastAPI, Biopython | 8000 | Track A step 1 — validate & translate DNA |
| `protein-analysis-service` | FastAPI, Biopython | 8001 | Track A step 2 — BLAST/UniProt/fold/structure |
| `drug-discovery-service` | FastAPI, PyTorch-Geometric, RDKit | 8002 | Shared output — pockets, GNN affinity, ranking |
| `histopathology-service` | FastAPI, torchvision | 8003 | Track B — WSI → mutation prediction + Grad-CAM |
| `postgres` | Postgres 16 | 5432 | Users, projects, jobs, reports |
| `redis` | Redis | 6379 | Per-service result caches |

A shared Docker volume `bio-data` is mounted at `/data/structures` in the
protein, drug, histopathology, and gateway containers. It carries PDB files,
Grad-CAM overlays, drug-docking visualizations, and uploaded WSI images between
services by path reference.

## Request lifecycle (a Job)

```
POST /api/v1/jobs                      (gateway, returns immediately, status=pending)
   └─ BackgroundTask: run_pipeline(job_id)
        Track A (if dna_sequence|gene_id):
          1. POST genomics  /api/v1/validate     -> validation_result
          2. pick longest ORF
          3. POST protein   /api/v1/analyze      -> analysis_result (incl. structure_3d)
          4. POST drug      /api/v1/predict      -> drug_result
        Track B (if wsi_image_path):
          5. POST histopath /api/v1/predict-mutation -> histopathology_result
        Fusion:
          6. run_fusion(analysis_result, histopathology_result) -> fusion_result
        status=completed
```

Each stage is wrapped in try/except: a single failing stage is recorded on the
job's `error` field and the pipeline continues, so a partial result is still
returned. The frontend polls `GET /api/v1/jobs/{id}` and renders progress from
`current_step`.

## Interface contracts

The canonical request/response shapes live in `packages/contracts/gtt_contracts`:

- `genomics.py` — `ValidationRequest` / `ValidationResponse` / `ORFRecord`
- `protein.py` — `AnalyzeRequest` / `AnalysisResponse` (BLAST, annotation, fold, `Structure3D`)
- `drug.py` — `PredictRequest` / `PredictResponse` (`Pocket`, `DrugCandidate`)
- `histopathology.py` — `HistopathologyRequest` / `HistopathologyResponse` (`MutationPrediction`)
- `fusion.py` — `FusionResult` / `GeneFusion`

Each microservice currently keeps a local copy of its models for historical
reasons; `gtt_contracts` is the source of truth and should be adopted by the
services over time to remove drift (e.g. `Structure3D` is defined slightly
differently in the protein vs. drug service today).

## Data model (gateway / Postgres)

```
User 1───* Project 1───* Job 1───* Report
User 1──────────────────* Job
```

- **User** — email, hashed password (bcrypt), role (`researcher` | `admin`), active flag.
- **Project** — grouping for analyses; has a cancer type.
- **Job** — one pipeline run: inputs + per-stage JSON results + lifecycle status.
- **Report** — a frozen summary snapshot of a completed job, renderable to HTML/PDF.

## Authentication

JWT bearer tokens (HS256). `POST /auth/register` (first user becomes admin),
`POST /auth/login`, and an OAuth2 `POST /auth/token` for Swagger. Protected
routes resolve the user via the `Authorization: Bearer` header. Admin-only
routes (user management) require `role == "admin"`.

## Known gaps / next steps

1. **Track B model weights** — the CNN+MIL run untrained; train ResNet50/ABMIL on
   TCGA WSIs and ship weights (mirroring the drug-discovery GNN approach).
2. **True genomic variant calling** in Track A would make the fusion layer
   clinically meaningful (today the genomic signal is a documented heuristic).
3. **Live 3D viewers** — the drug service emits py3Dmol HTML and the histo service
   emits Grad-CAM PNGs to the shared volume; the frontend currently references
   them by path. Serve them (static route or object storage) for in-browser viewing.
4. **DB migrations** — the gateway uses `create_all` on startup; add Alembic for
   schema evolution.
5. **Async job queue** — pipeline runs in a FastAPI BackgroundTask; move to Celery/RQ
   for retries, concurrency limits, and durability.
```
