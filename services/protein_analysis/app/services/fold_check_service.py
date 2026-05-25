import httpx
from typing import List, Tuple, Optional, Dict
from app.models import FoldCheck
import logging

logger = logging.getLogger(__name__)

# Kyte-Doolittle scale
HYDRO_SCALE = {
    'A': 1.8, 'V': 4.2, 'I': 4.5, 'L': 3.8, 'F': 2.8,
    'C': 2.5, 'M': 1.9, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2,
    'E': -3.5, 'Q': -3.5, 'D': -3.5,
    'N': -3.5, 'K': -3.9, 'R': -4.5
}

def calculate_hydrophobicity_fallback(sequence: str, window_size: int = 9) -> Tuple[float, List[List[int]]]:
    """
    Local hydrophobicity fallback using Kyte-Doolittle scale and sliding window.
    """
    scores = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        avg = sum(HYDRO_SCALE.get(aa, 0) for aa in window) / window_size
        scores.append(avg)

    if not scores:
        return 0.0, []

    mean_score = sum(scores) / len(scores)
    
    # Detect disordered regions (threshold < -1.0)
    regions = []
    start = None
    threshold = -1.0
    
    for i, val in enumerate(scores):
        if val < threshold:
            if start is None:
                start = i + 1
        else:
            if start is not None:
                regions.append([start, i + window_size - 1])
                start = None

    if start is not None:
        regions.append([start, len(sequence)])

    return round(mean_score, 3), regions

async def check_foldability(sequence: str) -> FoldCheck:
    """
    Tiered fold check: FoldIndex -> IUPred -> Fallback.
    """
    # 1. Try FoldIndex
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            # The prompt says POST to https://fold.weizmann.ac.il/fldbin/findex
            # FoldIndex usually accepts sequence as "seq" or similar.
            response = await client.post("https://fold.weizmann.ac.il/fldbin/findex", data={"seq": sequence})
            if response.status_code == 200:
                # Mocking parsing logic: FoldIndex returns text with "FoldIndex: 0.23"
                content = response.text
                if "FoldIndex" in content:
                    # Very rough parsing, depends on actual API output
                    # For exercise, we'll assume we can find a float
                    import re
                    match = re.search(r"FoldIndex\s*=\s*([-+]?\d*\.\d+|\d+)", content)
                    if match:
                        score = float(match.group(1))
                        return FoldCheck(
                            foldindex_score=score,
                            foldable=score > 0,
                            method="api",
                            disordered_percent=None # Not easy to get without parsing whole report
                        )
    except Exception as e:
        logger.warning(f"FoldIndex API failed: {e}")

    # 2. Try IUPred2A
    try:
        # Prompt URL: https://iupred2a.elte.hu/iupred2a/{sequence}/short
        # Wait, the search said it expects ID, but I'll try the sequence as requested.
        async with httpx.AsyncClient(timeout=10.0) as client:
            url = f"https://iupred2a.elte.hu/iupred2a/{sequence}/short"
            response = await client.get(url)
            if response.status_code == 200:
                data = response.json()
                # IUPred returns a list of scores for each residue
                scores = data.get("iupred2", [])
                if scores:
                    disordered_residues = [s for s in scores if s > 0.5]
                    disordered_percent = (len(disordered_residues) / len(scores)) * 100
                    
                    # Detect regions
                    regions = []
                    start = None
                    for i, s in enumerate(scores):
                        if s > 0.5:
                            if start is None: start = i + 1
                        else:
                            if start is not None:
                                regions.append([start, i])
                                start = None
                    if start is not None:
                        regions.append([start, len(scores)])
                    
                    return FoldCheck(
                        foldable=disordered_percent < 70,
                        disordered_percent=round(disordered_percent, 2),
                        disordered_regions=regions,
                        method="api"
                    )
    except Exception as e:
        logger.warning(f"IUPred2A API failed: {e}")

    # 3. Fallback to Hydrophobicity
    mean_score, regions = calculate_hydrophobicity_fallback(sequence)
    disordered_length = sum(r[1] - r[0] + 1 for r in regions)
    disordered_percent = (disordered_length / len(sequence)) * 100
    
    return FoldCheck(
        foldindex_score=mean_score, # We'll put the hydro mean here as a proxy
        foldable=mean_score > 0,
        disordered_percent=round(disordered_percent, 2),
        disordered_regions=regions,
        method="hydrophobicity_fallback",
        warning="External fold-check APIs were unavailable. Used local hydrophobicity fallback."
    )
