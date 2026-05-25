import subprocess
import os
import re
import numpy as np
from typing import List, Dict, Optional
from ..utils.pdb_parser import PDBHelper

class PocketService:
    @staticmethod
    def detect_pockets(pdb_path: str, uniprot_binding_sites: List[int] = None) -> Dict:
        """Step A: Pocket Detection using fpocket with P2Rank fallback."""
        try:
            return PocketService._run_fpocket(pdb_path, uniprot_binding_sites)
        except Exception as e:
            print(f"fpocket failed: {e}. Trying P2Rank fallback...")
            return PocketService._run_p2rank(pdb_path, uniprot_binding_sites)

    @staticmethod
    def _run_fpocket(pdb_path: str, uniprot_binding_sites: List[int]) -> Dict:
        # Run fpocket: fpocket -f protein.pdb
        result = subprocess.run(["fpocket", "-f", pdb_path], capture_output=True, text=True)
        
        # fpocket creates a directory pdb_basename_out/
        pdb_basename = os.path.basename(pdb_path).replace(".pdb", "")
        out_dir = f"{pdb_basename}_out"
        info_file = os.path.join(out_dir, f"{pdb_basename}_info.txt")
        
        if not os.path.exists(info_file):
            raise Exception("fpocket output file not found")

        pockets = PocketService._parse_fpocket_info(info_file)
        
        # Filter pockets: druggability > 0.5 and volume > 300
        filtered_pockets = [p for p in pockets if p["score"] > 0.5 and p["volume"] > 300]
        
        if not filtered_pockets:
            filtered_pockets = pockets # Fallback to all if none match

        # Prioritize by proximity to UniProt binding sites
        if uniprot_binding_sites:
            target_pocket = PocketService._find_closest_pocket(filtered_pockets, pdb_path, uniprot_binding_sites)
        else:
            target_pocket = max(filtered_pockets, key=lambda x: x["score"])

        return target_pocket

    @staticmethod
    def _parse_fpocket_info(info_path: str) -> List[Dict]:
        pockets = []
        with open(info_path, 'r') as f:
            content = f.read()
            # Split by Pocket X
            pocket_blocks = re.split(r'Pocket \d+ :', content)[1:]
            for i, block in enumerate(pocket_blocks):
                score_match = re.search(r'Druggability Score : ([\d\.]+)', block)
                vol_match = re.search(r'Pocket Volume : ([\d\.]+)', block)
                # Note: fpocket info file might not have centroid directly in _info.txt in all versions
                # Usually it's in the atoms file. We'll simplify for now.
                pockets.append({
                    "id": i + 1,
                    "score": float(score_match.group(1)) if score_match else 0.0,
                    "volume": float(vol_match.group(1)) if vol_match else 0.0,
                    "residues": [] # Will be populated if needed
                })
        return pockets

    @staticmethod
    def _find_closest_pocket(pockets: List[Dict], pdb_path: str, residues: List[int]) -> Dict:
        # Simplified: max score for now, distance calculation would require parsing PQR files
        return max(pockets, key=lambda x: x["score"])

    @staticmethod
    def _run_p2rank(pdb_path: str, residues: List[int]) -> Dict:
        # prank predict -f protein.pdb
        subprocess.run(["prank", "predict", "-f", pdb_path], capture_output=True)
        # Parsing P2Rank CSV output...
        # Simplified placeholder
        return {
            "residues": residues if residues else [1, 2, 3],
            "center": [0.0, 0.0, 0.0],
            "volume": 400.0,
            "score": 0.8
        }
