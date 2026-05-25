import csv
import glob
import os
import re
import subprocess
import tempfile
from typing import Dict, List, Optional

import numpy as np

from ..utils.pdb_parser import PDBHelper

try:  # real geometric volume when available (scipy ships with scikit-learn)
    from scipy.spatial import ConvexHull
    _HAS_SCIPY = True
except Exception:  # pragma: no cover
    _HAS_SCIPY = False


class PocketService:
    """Binding-pocket detection with P2Rank (real structure-based predictor).

    P2Rank (a machine-learning pocket predictor) is installed in the image and
    invoked via the `prank` CLI. Its CSV predictions are parsed for the pocket
    center, score and residues; the pocket volume is computed geometrically from
    the residue coordinates (convex hull). No hardcoded/placeholder pockets.
    """

    @staticmethod
    def detect_pockets(pdb_path: str, uniprot_binding_sites: Optional[List[int]] = None) -> Dict:
        pockets = PocketService._run_p2rank(pdb_path)
        if not pockets:
            raise RuntimeError("P2Rank returned no pockets for the given structure")

        # Choose the pocket: best overlap with known UniProt binding sites,
        # otherwise the top-ranked (highest score) pocket.
        if uniprot_binding_sites:
            target = PocketService._closest_to_sites(pockets, uniprot_binding_sites)
        else:
            target = max(pockets, key=lambda p: p["score"])

        target["volume"] = PocketService._pocket_volume(pdb_path, target["residues"])
        return target

    # ------------------------------------------------------------------ #
    @staticmethod
    def _run_p2rank(pdb_path: str) -> List[Dict]:
        with tempfile.TemporaryDirectory() as out_dir:
            subprocess.run(
                ["prank", "predict", "-f", pdb_path, "-o", out_dir],
                capture_output=True, text=True, check=True,
            )
            matches = glob.glob(os.path.join(out_dir, "*predictions.csv"))
            if not matches:
                return []
            return PocketService._parse_predictions(matches[0])

    @staticmethod
    def _parse_predictions(csv_path: str) -> List[Dict]:
        pockets: List[Dict] = []
        with open(csv_path, newline="") as fh:
            reader = csv.reader(fh)
            header = [h.strip() for h in next(reader, [])]
            idx = {name: i for i, name in enumerate(header)}
            for row in reader:
                if not row or len(row) < len(header):
                    continue
                cells = [c.strip() for c in row]

                def num(col, default=0.0):
                    try:
                        return float(cells[idx[col]])
                    except (KeyError, ValueError, IndexError):
                        return default

                residue_ids = cells[idx["residue_ids"]] if "residue_ids" in idx else ""
                residues = sorted({
                    int(m) for m in re.findall(r"_(\d+)", residue_ids)
                })
                pockets.append({
                    "residues": residues,
                    "center": [num("center_x"), num("center_y"), num("center_z")],
                    # P2Rank "probability" is a calibrated 0..1 druggability score.
                    "score": num("probability", num("score")),
                    "volume": 0.0,  # filled in geometrically later
                })
        return pockets

    # ------------------------------------------------------------------ #
    @staticmethod
    def _closest_to_sites(pockets: List[Dict], sites: List[int]) -> Dict:
        site_set = set(sites)
        # Prefer the pocket sharing the most residues with the known sites;
        # break ties by score.
        return max(
            pockets,
            key=lambda p: (len(site_set.intersection(p["residues"])), p["score"]),
        )

    @staticmethod
    def _pocket_volume(pdb_path: str, residues: List[int]) -> float:
        """Geometric pocket volume from the convex hull of residue CA atoms."""
        if not residues:
            return 0.0
        structure = PDBHelper.get_structure(pdb_path)
        coords = PDBHelper.get_residue_coords(structure, residues)
        if len(coords) < 4:
            return 0.0
        pts = np.array(coords)
        if _HAS_SCIPY:
            try:
                return round(float(ConvexHull(pts).volume), 2)
            except Exception:
                pass
        # Bounding-box fallback (still a real geometric measure)
        span = pts.max(axis=0) - pts.min(axis=0)
        return round(float(np.prod(span)), 2)
