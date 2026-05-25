# gtt-contracts

Canonical interface contracts for the **Genomics-to-Therapy AI Platform**.

This package is the single source of truth for the request/response shapes that
travel between microservices and the API gateway. Each microservice currently
keeps its own local Pydantic models (historical reasons); the gateway imports
**these** models when it orchestrates the pipeline, so they define the contract.

| Module | Service | Endpoint |
|--------|---------|----------|
| `genomics.py` | Genomics Validation (Track A) | `POST /api/v1/validate` |
| `protein.py` | Protein Analysis (Track A) | `POST /api/v1/analyze` |
| `drug.py` | Drug Discovery (shared output) | `POST /api/v1/predict` |
| `histopathology.py` | Histopathology (Track B) | `POST /api/v1/predict-mutation` |
| `fusion.py` | Gateway fusion layer | cross-modal mutation validation |

Install (editable) for local dev:

```bash
pip install -e packages/contracts
```

Then `from gtt_contracts.protein import AnalysisResponse`.
