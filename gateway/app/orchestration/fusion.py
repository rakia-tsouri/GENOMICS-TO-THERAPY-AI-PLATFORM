"""Cross-modal fusion layer.

Reconciles mutation evidence from Track A (genomic sequence analysis) and
Track B (histopathology visual prediction). Concordant signals raise
confidence; discordant signals are flagged for manual review.

NOTE: Track A as currently built does sequence validation + protein/structure
analysis, not variant calling. The genomic signal here is therefore a
documented heuristic ("is this gene the one analyzed, and is it disease-
associated?"). Wiring true variant calling into Track A would make this
fusion clinically meaningful — the structure below is ready for that.
"""
from __future__ import annotations

from typing import Optional


def _genomic_signal_for_gene(gene: str, analysis_result: Optional[dict]) -> Optional[bool]:
    """Heuristic Track-A signal for a gene. Returns None if not assessed."""
    if not analysis_result:
        return None
    analyzed_gene = (analysis_result.get("gene_id") or "").upper()
    if gene.upper() not in analyzed_gene:
        return None  # Track A did not analyze this gene
    annotation = analysis_result.get("annotation") or {}
    diseases = annotation.get("diseases") or []
    # Disease-associated annotation is treated as genomic concern for the gene.
    return bool(diseases)


def _confidence_label_to_weight(label: str) -> float:
    return {"high": 0.9, "medium": 0.6, "low": 0.3}.get(label, 0.5)


def run_fusion(
    target_mutations: list[str],
    analysis_result: Optional[dict],
    histopathology_result: Optional[dict],
) -> dict:
    """Return a FusionResult-shaped dict (see gtt_contracts.fusion)."""
    visual_by_gene: dict[str, dict] = {}
    if histopathology_result:
        for m in histopathology_result.get("mutations", []):
            visual_by_gene[m["gene"].upper()] = m

    genes_out = []
    flagged = []
    both_signal_pairs = 0
    agreeing_pairs = 0

    for gene in target_mutations:
        gkey = gene.upper()
        visual = visual_by_gene.get(gkey)
        visual_signal = visual["mutated"] if visual else None
        visual_prob = visual["probability"] if visual else None
        genomic_signal = _genomic_signal_for_gene(gene, analysis_result)

        if genomic_signal is not None and visual_signal is not None:
            both_signal_pairs += 1
            if genomic_signal == visual_signal:
                agreement = "concordant"
                agreeing_pairs += 1
                # Reinforce: average of modality confidences, boosted.
                base = (
                    _confidence_label_to_weight(visual.get("confidence", "medium"))
                    + (0.8 if genomic_signal else 0.5)
                ) / 2
                combined = min(1.0, base + 0.15)
                flag = False
                note = "Both modalities agree."
            else:
                agreement = "discordant"
                combined = 0.35
                flag = True
                note = "Genomic and visual evidence disagree — manual review recommended."
                flagged.append(gene)
        elif visual_signal is not None:
            agreement = "single_modality"
            combined = _confidence_label_to_weight(visual.get("confidence", "medium"))
            flag = False
            note = "Visual (Track B) evidence only."
        elif genomic_signal is not None:
            agreement = "single_modality"
            combined = 0.7 if genomic_signal else 0.4
            flag = False
            note = "Genomic (Track A) evidence only."
        else:
            agreement = "single_modality"
            combined = 0.0
            flag = False
            note = "No evidence from either modality."

        genes_out.append({
            "gene": gene,
            "genomic_signal": genomic_signal,
            "visual_signal": visual_signal,
            "visual_probability": visual_prob,
            "agreement": agreement,
            "combined_confidence": round(combined, 3),
            "flag_for_review": flag,
            "note": note,
        })

    # Overall agreement = worst case across genes.
    if any(g["agreement"] == "discordant" for g in genes_out):
        overall = "discordant"
    elif any(g["agreement"] == "concordant" for g in genes_out):
        overall = "concordant"
    else:
        overall = "single_modality"

    # Simple observed-agreement proxy for Cohen's kappa across dual-signal genes.
    kappa = round(agreeing_pairs / both_signal_pairs, 3) if both_signal_pairs else None

    return {
        "genes": genes_out,
        "overall_agreement": overall,
        "cohen_kappa": kappa,
        "flagged_genes": flagged,
    }
