"""Thin async HTTP clients for the downstream microservices.

Each function calls one service endpoint and returns the parsed JSON dict.
Network/HTTP errors raise httpx exceptions; the orchestrator decides how to
degrade gracefully.
"""
import httpx

from ..config import settings


async def _post(url: str, payload: dict) -> dict:
    async with httpx.AsyncClient(timeout=settings.service_timeout_sec) as client:
        resp = await client.post(url, json=payload)
        resp.raise_for_status()
        return resp.json()


async def call_validate(payload: dict) -> dict:
    """Genomics Validation Service — POST /api/v1/validate"""
    return await _post(f"{settings.genomics_url}/api/v1/validate", payload)


async def call_analyze(payload: dict) -> dict:
    """Protein Analysis Service — POST /api/v1/analyze"""
    return await _post(f"{settings.protein_url}/api/v1/analyze", payload)


async def call_predict_drugs(payload: dict) -> dict:
    """Drug Discovery Service — POST /api/v1/predict"""
    return await _post(f"{settings.drug_url}/api/v1/predict", payload)


async def call_predict_mutation(payload: dict) -> dict:
    """Histopathology Service (Track B) — POST /api/v1/predict-mutation"""
    return await _post(f"{settings.histopathology_url}/api/v1/predict-mutation", payload)
