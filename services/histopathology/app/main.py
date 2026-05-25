"""Histopathology Visual Analysis Service (Track B).

POST /api/v1/predict-mutation : WSI -> patches -> CNN features -> ABMIL ->
per-gene mutation prediction + Grad-CAM overlay.
GET  /health

Pipeline stages are each wrapped in try/except (mirroring
protein_analysis/app/main.py) so a single stage failure degrades the request
gracefully (status / reason / warnings) instead of returning a 500.

The model runs on UNTRAINED weights; a fixed warning is always appended so
downstream consumers know predictions are illustrative scaffolding.
"""
import time
import logging

from fastapi import FastAPI

from app.models import HistopathologyRequest, HistopathologyResponse, MutationPrediction
from app.services import patch_service
from app.services.cnn_service import cnn_feature_extractor
from app.services.mil_service import mil_service, confidence_from_prob

# Configure logging (plain logging, same style as the other services).
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Histopathology Service",
    description=(
        "Track B microservice: predict TP53/IDH1/KRAS mutations directly from "
        "WSI tissue morphology (ResNet50 + attention-MIL) with Grad-CAM "
        "explainability."
    ),
    version="1.0.0",
)

# Always-on disclaimer: there are no trained CNN/MIL weights.
UNTRAINED_WARNING = (
    "Histopathology model is running on untrained weights; predictions are "
    "illustrative scaffolding, not validated."
)


@app.post("/api/v1/predict-mutation", response_model=HistopathologyResponse)
async def predict_mutation(input_data: HistopathologyRequest):
    """Full histopathology mutation-prediction pipeline (error-tolerant)."""
    start_time = time.time()
    warnings = [UNTRAINED_WARNING]
    logger.info(
        "Starting histopathology analysis for WSI: %s (genes=%s)",
        input_data.wsi_image_path,
        input_data.target_mutations,
    )

    num_patches = 0
    patches_kept = 0
    patches = []

    # ----- Stage A: patch extraction (Otsu + tissue filter + stain norm) -----
    try:
        patch_result = patch_service.extract_patches(
            input_data.wsi_image_path,
            patch_size=input_data.patch_size,
            magnification=input_data.magnification,
        )
        patches = patch_result["patches"]
        num_patches = patch_result["num_patches"]
        patches_kept = patch_result["patches_kept"]
    except Exception as e:
        logger.error("Error in patch extraction stage: %s", e)
        # Unrecoverable: without patches the whole pipeline cannot proceed.
        return HistopathologyResponse(
            status="failed",
            reason=f"patch_extraction_failed: {e}",
            num_patches=num_patches,
            patches_kept=patches_kept,
            warnings=warnings,
            processing_time_sec=round(time.time() - start_time, 2),
        )

    if patches_kept == 0:
        logger.warning("No tissue patches kept after filtering for %s", input_data.wsi_image_path)
        warnings.append("No tissue patches passed the Otsu/background filter.")
        return HistopathologyResponse(
            status="failed",
            reason="no_tissue_patches",
            num_patches=num_patches,
            patches_kept=patches_kept,
            warnings=warnings,
            processing_time_sec=round(time.time() - start_time, 2),
        )

    # ----- Stage B: CNN feature extraction (ResNet50) -----
    features = None
    try:
        features = cnn_feature_extractor.extract_features(patches)
    except Exception as e:
        logger.error("Error in CNN feature extraction stage: %s", e)
        warnings.append(f"cnn_feature_extraction_failed: {e}")

    # ----- Stage C: ABMIL aggregation + per-gene prediction -----
    mutations = []
    attention = None
    if features is not None and len(features) > 0:
        try:
            mil_out = mil_service.predict(features, input_data.target_mutations)
            attention = mil_out["attention"]
            for p in mil_out["predictions"]:
                mutations.append(
                    MutationPrediction(
                        gene=p["gene"],
                        mutated=p["mutated"],
                        probability=p["probability"],
                        confidence=p["confidence"],
                    )
                )
        except Exception as e:
            logger.error("Error in MIL prediction stage: %s", e)
            warnings.append(f"mil_prediction_failed: {e}")
    else:
        warnings.append("Skipping MIL stage: no CNN features available.")

    # If MIL produced nothing, emit safe neutral predictions so the response
    # still contains one MutationPrediction per requested gene.
    if not mutations:
        for gene in input_data.target_mutations:
            mutations.append(
                MutationPrediction(
                    gene=gene,
                    mutated=False,
                    probability=0.5,
                    confidence=confidence_from_prob(0.5),
                )
            )

    # ----- Stage D: Grad-CAM overlay for the highest-attention patch -----
    gradcam_path = None
    try:
        # Lazy import keeps cv2/torch import cost out of the module-load path
        # and isolates any cv2 graphics-lib issues to this stage.
        from app.services import gradcam_service

        gradcam_path = gradcam_service.generate_overlay(
            patches=patches,
            attention=attention,
            cnn_extractor=cnn_feature_extractor,
        )
    except Exception as e:
        logger.error("Error in Grad-CAM stage: %s", e)
        warnings.append(f"gradcam_failed: {e}")

    processing_time = time.time() - start_time
    logger.info("Histopathology analysis completed in %.2f s.", processing_time)

    # Status is "success" as long as we return per-gene predictions, even if the
    # CNN/MIL/Grad-CAM stages degraded (recorded in warnings).
    return HistopathologyResponse(
        status="success",
        num_patches=num_patches,
        patches_kept=patches_kept,
        mutations=mutations,
        gradcam_overlay_path=gradcam_path,
        model="ResNet50+ABMIL",
        processing_time_sec=round(processing_time, 2),
        warnings=warnings,
    )


@app.get("/health")
async def health_check():
    return {"status": "healthy"}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8003)
