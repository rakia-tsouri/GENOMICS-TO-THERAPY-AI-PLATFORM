# Genomics-to-Therapy AI Platform

A unified, multi-modal biomedical AI pipeline that goes **from a patient's DNA
sequence and/or tumor tissue image to ranked drug candidates** — cross-validating
mutations across genomics and histopathology along the way.

It merges two AI tracks behind one product:

- **Track A — Genomic Sequence Analysis:** DNA → validation → translation →
  BLAST/UniProt annotation → foldability → 3D protein structure.
- **Track B — Histopathology Visual Analysis:** WSI tissue image → patch
  extraction → CNN features → MIL aggregation → mutation prediction → Grad-CAM.
- **Fusion Layer:** reconciles mutation evidence from both tracks (concordant →
  higher confidence; discordant → flagged for review).
- **Shared output — Drug Discovery:** binding-pocket detection → molecular
  graphs → GNN affinity prediction → ranked candidates + report.

```
                    ┌────────────────────────────── Frontend (Next.js, :3000) ──────────────────────────────┐
                    │   auth · dashboard · projects · analyses · results · reports · admin users             │
                    └───────────────────────────────────────────┬───────────────────────────────────────────┘
                                                                 │  REST + JWT
                                          ┌──────────────────────▼───────────────────────┐
                                          │      API Gateway (FastAPI, :8080)             │
                                          │  auth · users · projects · jobs · reports     │
                                          │  orchestration · cross-modal fusion           │
                                          └───┬───────────┬───────────┬───────────┬───────┘
                       Track A               │           │           │           │      Track B
        ┌─────────────────────────┬─────────▼──┐   ┌─────▼──────┐  ┌─▼──────────┐ │ ┌────────────────────┐
        │ genomics-validation :8000│  protein-  │   │  drug-     │  │            │ └►│ histopathology :8003│
        │  /validate               │  analysis  │   │ discovery  │  │            │   │ /predict-mutation   │
        │                          │  :8001     │   │  :8002     │  │            │   │ (WSI→CNN→MIL→GradCAM)│
        └──────────────────────────┴───/analyze─┘   └──/predict──┘  └────────────┘   └────────────────────┘
                                          │                                  Postgres (:5432) · Redis (:6379)
```

## Repository layout

```
.
├── docker-compose.yml          # full stack: db, cache, 4 services, gateway, frontend
├── .env.example                # copy to .env
├── docs/architecture.md        # deeper architecture notes
├── packages/
│   └── contracts/              # gtt-contracts: canonical Pydantic interface contracts
├── services/
│   ├── genomics_validation/    # Track A · POST /api/v1/validate            (8000)
│   ├── protein_analysis/       # Track A · POST /api/v1/analyze             (8001)
│   ├── drug_discovery/         # shared  · POST /api/v1/predict + GNN train (8002)
│   └── histopathology/         # Track B · POST /api/v1/predict-mutation    (8003)
├── gateway/                    # FastAPI gateway: auth + Postgres + orchestration + fusion (8080)
└── frontend/                   # Next.js + TS + Tailwind product UI         (3000)
```

## Service status

| Component | Status | Notes |
|-----------|--------|-------|
| Genomics Validation (Track A) | ✅ **Implemented** | FASTA/raw/Gene-ID input, NCBI fetch, ORF, translation, Redis cache |
| Protein Analysis (Track A) | ✅ **Implemented** | BLAST, UniProt, FoldIndex/IUPred, AlphaFold/PDB/ESMFold |
| Drug Discovery (shared) | ✅ **Implemented** | Pockets, ChEMBL/ZINC, RDKit graphs, **trained GNN**, ADMET, viz |
| Histopathology (Track B) | 🟡 **Scaffolded** | Real pipeline structure; CNN/MIL run on **untrained weights** (clearly flagged). Needs trained ResNet/ABMIL weights on TCGA. |
| API Gateway | 🟡 **New (this pass)** | Auth, users, projects, jobs, reports, orchestration, fusion. Needs end-to-end run. |
| Frontend | 🟡 **New (this pass)** | Next.js app with full UX. Run `npm install` first. |
| Fusion layer | 🟡 **Heuristic** | Works; production cross-validation needs variant calling wired into Track A. |

## Quick start (Docker)

```bash
cp .env.example .env          # adjust secrets
docker compose up --build
```

Then open:

- Frontend: http://localhost:3000
- Gateway API docs (Swagger): http://localhost:8080/docs
- A default admin is seeded from `.env` (`ADMIN_EMAIL` / `ADMIN_PASSWORD`).

> **Heads up:** the `drug-discovery` and `histopathology` images pull large ML
> deps (torch, torch-geometric, rdkit, torchvision) — the first build is slow.
> The bioinformatics services also call external APIs (NCBI, BLAST, UniProt,
> AlphaFold, ESMFold, ChEMBL) at runtime, so they need internet access.

## Local development (without Docker)

Each service is an independent FastAPI app:

```bash
# Gateway (boots on SQLite by default — no Postgres needed for quick dev)
cd gateway && pip install -r requirements.txt && uvicorn app.main:app --reload --port 8080

# A Track A service, e.g. genomics
cd services/genomics_validation && pip install -r requirements.txt
PYTHONPATH=. uvicorn api.main:app --reload --port 8000

# Frontend
cd frontend && npm install && npm run dev
```

The shared interface contracts live in `packages/contracts` (`pip install -e packages/contracts`).

## How a job flows

1. User signs in, creates a **Project**, then starts an **Analysis (Job)** with a
   DNA sequence / Gene ID and/or an uploaded WSI image.
2. The gateway runs the pipeline in the background, updating `status`/`current_step`:
   `validating → analyzing → predicting_drugs → histopathology → fusion → done`.
3. The results page polls the job and renders each stage; the user can generate a
   **Report** (viewable in-app and exportable to HTML/PDF).

## Credits

Track A microservices (genomics, protein analysis, drug discovery incl. the GNN
training pipeline) by **rakia-tsouri**. Monorepo reorganization, API gateway,
Track B scaffold, and frontend added in a later pass.

> For research and educational use only. Not a medical device.
