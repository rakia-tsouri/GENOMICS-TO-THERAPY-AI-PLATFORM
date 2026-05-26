"""Build a structured report summary from a completed job and render it to HTML.

The HTML is intentionally self-contained (inline CSS) so it can be downloaded
and opened anywhere, or piped to a PDF engine later.
"""
from __future__ import annotations

import html
from typing import Any

from ..models import Job


def build_summary(job: Job) -> dict[str, Any]:
    """Flatten the job's stage results into a compact, render-ready dict."""
    analysis = job.analysis_result or {}
    drug = job.drug_result or {}
    histo = job.histopathology_result or {}
    fusion = job.fusion_result or {}

    top_candidates = drug.get("top_candidates", [])[:5]

    return {
        "job_id": job.id,
        "job_name": job.name,
        "created_at": job.created_at.isoformat() if job.created_at else None,
        "inputs": {
            "gene_id": job.gene_id,
            "has_dna": bool(job.dna_sequence),
            "has_wsi": bool(job.wsi_image_path),
            "target_mutations": job.target_mutations,
        },
        "genomics": {
            "valid": (job.validation_result or {}).get("valid"),
            "gc_percent": (job.validation_result or {}).get("gc_percent"),
            "orf_count": len((job.validation_result or {}).get("orfs", [])),
        },
        "protein": {
            "blast_status": (analysis.get("blast") or {}).get("protein_status"),
            "top_hit": (analysis.get("blast") or {}).get("top_hit_name"),
            "uniprot_id": (analysis.get("annotation") or {}).get("uniprot_id"),
            "foldable": (analysis.get("fold_check") or {}).get("foldable"),
            "structure_source": (analysis.get("structure_3d") or {}).get("source"),
            "plddt_mean": (analysis.get("structure_3d") or {}).get("plddt_mean"),
        },
        "drugs": {
            "pocket_score": (drug.get("pocket") or {}).get("score"),
            "candidate_count": len(drug.get("top_candidates", [])),
            "top_candidates": [
                {
                    "chembl_id": c.get("chembl_id"),
                    "smiles": c.get("smiles"),
                    "binding_score": c.get("binding_score"),
                    "toxicity_risk": c.get("toxicity_risk"),
                }
                for c in top_candidates
            ],
        },
        "histopathology": {
            "model": histo.get("model"),
            "mutations": histo.get("mutations", []),
            "patches_kept": histo.get("patches_kept"),
        },
        "fusion": fusion,
        "warnings": list(
            {
                *(job.validation_result or {}).get("warnings", []),
                *((analysis.get("fold_check") or {}).get("warning") and
                  [(analysis.get("fold_check") or {}).get("warning")] or []),
                *histo.get("warnings", []),
            }
        ),
    }


def _row(label: str, value: Any) -> str:
    return f"<tr><th>{html.escape(str(label))}</th><td>{html.escape(str(value))}</td></tr>"


def render_html(summary: dict[str, Any]) -> str:
    g, p, d = summary["genomics"], summary["protein"], summary["drugs"]
    fusion = summary.get("fusion") or {}

    drug_rows = "".join(
        f"<tr><td>{html.escape(str(c.get('chembl_id') or '—'))}</td>"
        f"<td class='mono'>{html.escape(str(c.get('smiles') or ''))}</td>"
        f"<td>{c.get('binding_score')}</td>"
        f"<td>{html.escape(str(c.get('toxicity_risk') or ''))}</td></tr>"
        for c in d.get("top_candidates", [])
    ) or "<tr><td colspan='4'>No drug candidates.</td></tr>"

    fusion_rows = "".join(
        f"<tr><td>{html.escape(str(gene.get('gene')))}</td>"
        f"<td>{gene.get('genomic_signal')}</td>"
        f"<td>{gene.get('visual_signal')}</td>"
        f"<td>{html.escape(str(gene.get('agreement')))}</td>"
        f"<td>{gene.get('combined_confidence')}</td>"
        f"<td>{'⚠️' if gene.get('flag_for_review') else ''}</td></tr>"
        for gene in fusion.get("genes", [])
    ) or "<tr><td colspan='6'>No fusion data.</td></tr>"

    warnings = "".join(f"<li>{html.escape(str(w))}</li>" for w in summary.get("warnings", []))

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Report — {html.escape(summary['job_name'])}</title>
<style>
 body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:40px;color:#1a1a2e;}}
 h1{{margin-bottom:0}} .sub{{color:#666;margin-top:4px}}
 h2{{margin-top:32px;border-bottom:2px solid #6c5ce7;padding-bottom:6px}}
 table{{border-collapse:collapse;width:100%;margin-top:8px}}
 th,td{{border:1px solid #ddd;padding:8px;text-align:left;font-size:14px;vertical-align:top}}
 th{{background:#f4f3ff;width:240px}}
 .mono{{font-family:ui-monospace,Menlo,monospace;font-size:12px;word-break:break-all}}
 .badge{{display:inline-block;padding:2px 10px;border-radius:12px;background:#6c5ce7;color:#fff;font-size:12px}}
 ul{{margin-top:8px}}
</style></head><body>
 <h1>Genomics-to-Therapy Report</h1>
 <p class="sub">{html.escape(summary['job_name'])} · Job #{summary['job_id']} · {html.escape(str(summary.get('created_at') or ''))}</p>
 <span class="badge">Gene: {html.escape(str(summary['inputs'].get('gene_id') or 'n/a'))}</span>

 <h2>Genomics (Track A)</h2>
 <table>{_row('Valid sequence', g.get('valid'))}{_row('GC %', g.get('gc_percent'))}{_row('ORFs found', g.get('orf_count'))}</table>

 <h2>Protein & Structure</h2>
 <table>
  {_row('BLAST status', p.get('blast_status'))}{_row('Top hit', p.get('top_hit'))}
  {_row('UniProt ID', p.get('uniprot_id'))}{_row('Foldable', p.get('foldable'))}
  {_row('Structure source', p.get('structure_source'))}{_row('Mean pLDDT', p.get('plddt_mean'))}
 </table>

 <h2>Drug Candidates</h2>
 <table><tr><th>ChEMBL ID</th><th>SMILES</th><th>Binding score</th><th>Toxicity</th></tr>{drug_rows}</table>

 <h2>Cross-Modal Fusion (Track A ⨯ Track B)</h2>
 <p>Overall agreement: <strong>{html.escape(str(fusion.get('overall_agreement', 'n/a')))}</strong>
    &nbsp; κ = {fusion.get('cohen_kappa')}</p>
 <table><tr><th>Gene</th><th>Genomic</th><th>Visual</th><th>Agreement</th><th>Confidence</th><th>Review</th></tr>{fusion_rows}</table>

 <h2>Warnings</h2>
 <ul>{warnings or '<li>None</li>'}</ul>

 <p class="sub" style="margin-top:40px">Generated by the Genomics-to-Therapy AI Platform. For research use only.</p>
</body></html>"""
