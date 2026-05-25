"""End-to-end pipeline orchestrator.

Runs as a background task after a Job is created. Chains:

  Track A:  validate -> analyze -> predict (drugs)
  Track B:  predict-mutation (histopathology)
  Fusion:   cross-modal reconciliation

Each stage is isolated: a failure in one stage is recorded on the job and the
pipeline continues with whatever evidence is available, so a partial result is
still useful (mirrors the error-tolerant style of the microservices).
"""
import logging
from datetime import datetime, timezone

from ..database import SessionLocal
from ..models import Job
from ..clients import services
from .fusion import run_fusion

logger = logging.getLogger("pipeline")


def _set(db, job: Job, **fields):
    for k, v in fields.items():
        setattr(job, k, v)
    db.add(job)
    db.commit()
    db.refresh(job)


def _pick_best_orf(validation_result: dict) -> dict | None:
    orfs = validation_result.get("orfs") or []
    return orfs[0] if orfs else None  # service returns ORFs sorted by length desc


async def run_pipeline(job_id: int) -> None:
    db = SessionLocal()
    try:
        job = db.query(Job).filter(Job.id == job_id).first()
        if not job:
            logger.error("Job %s vanished before pipeline start", job_id)
            return

        _set(db, job, status="running", current_step="starting", error="")

        analysis_result = None
        validation_result = None

        # ---------------- Track A: genomics ----------------
        if job.dna_sequence or job.gene_id:
            try:
                _set(db, job, current_step="validating")
                validation_result = await services.call_validate(
                    {"dna_sequence": job.dna_sequence, "gene_id": job.gene_id}
                )
                _set(db, job, validation_result=validation_result)
            except Exception as e:
                logger.exception("validate failed")
                _set(db, job, error=f"validation: {e}")

            best_orf = _pick_best_orf(validation_result) if validation_result else None
            if validation_result and validation_result.get("valid") and best_orf:
                try:
                    _set(db, job, current_step="analyzing")
                    analyze_payload = {
                        "valid": True,
                        "gene_id": validation_result.get("gene_id") or job.gene_id or "Unknown",
                        "dna_length": validation_result.get("dna_length", 0),
                        "protein_seq": best_orf["protein_seq"],
                        "protein_length": best_orf["protein_length"],
                        "gc_percent": validation_result.get("gc_percent", 0.0),
                        "warnings": validation_result.get("warnings", []),
                    }
                    analysis_result = await services.call_analyze(analyze_payload)
                    _set(db, job, analysis_result=analysis_result)
                except Exception as e:
                    logger.exception("analyze failed")
                    _set(db, job, error=f"{job.error} | analysis: {e}".strip(" |"))

            # ---------------- Drug discovery (shared output) ----------------
            structure = (analysis_result or {}).get("structure_3d") or {}
            if structure.get("pdb_file_path"):
                try:
                    _set(db, job, current_step="predicting_drugs")
                    annotation = (analysis_result or {}).get("annotation") or {}
                    predict_payload = {
                        "gene_id": analysis_result.get("gene_id"),
                        "protein_seq": analysis_result.get("protein_seq", ""),
                        "structure_3d": {
                            "pdb_file_path": structure["pdb_file_path"],
                            "plddt_mean": structure.get("plddt_mean") or 0.0,
                            "plddt_per_residue": structure.get("plddt_per_residue", []),
                        },
                        "annotation": {
                            "uniprot_id": annotation.get("uniprot_id"),
                            "binding_sites": annotation.get("binding_sites", []),
                        },
                    }
                    drug_result = await services.call_predict_drugs(predict_payload)
                    _set(db, job, drug_result=drug_result)
                except Exception as e:
                    logger.exception("drug prediction failed")
                    _set(db, job, error=f"{job.error} | drug: {e}".strip(" |"))

        # ---------------- Track B: histopathology ----------------
        histo_result = None
        if job.wsi_image_path:
            try:
                _set(db, job, current_step="histopathology")
                histo_result = await services.call_predict_mutation(
                    {
                        "wsi_image_path": job.wsi_image_path,
                        "cancer_type": None,
                        "target_mutations": job.target_mutations or ["TP53", "IDH1", "KRAS"],
                    }
                )
                _set(db, job, histopathology_result=histo_result)
            except Exception as e:
                logger.exception("histopathology failed")
                _set(db, job, error=f"{job.error} | histopathology: {e}".strip(" |"))

        # ---------------- Fusion ----------------
        try:
            _set(db, job, current_step="fusion")
            fusion = run_fusion(
                job.target_mutations or ["TP53", "IDH1", "KRAS"],
                analysis_result,
                histo_result,
            )
            _set(db, job, fusion_result=fusion)
        except Exception as e:
            logger.exception("fusion failed")
            _set(db, job, error=f"{job.error} | fusion: {e}".strip(" |"))

        _set(
            db,
            job,
            status="completed",
            current_step="done",
            completed_at=datetime.now(timezone.utc),
        )
        logger.info("Job %s completed", job_id)

    except Exception as e:  # pragma: no cover - safety net
        logger.exception("pipeline crashed")
        job = db.query(Job).filter(Job.id == job_id).first()
        if job:
            _set(db, job, status="failed", error=str(e))
    finally:
        db.close()
