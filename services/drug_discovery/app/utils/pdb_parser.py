from Bio.PDB import PDBParser
import numpy as np

class PDBHelper:
    @staticmethod
    def get_structure(pdb_path: str):
        parser = PDBParser(QUIET=True)
        return parser.get_structure("protein", pdb_path)

    @staticmethod
    def calculate_centroid(residues_coords):
        if not residues_coords:
            return [0.0, 0.0, 0.0]
        return np.mean(residues_coords, axis=0).tolist()

    @staticmethod
    def get_residue_coords(structure, residue_ids):
        coords = []
        for residue in structure.get_residues():
            if residue.get_id()[1] in residue_ids:
                if 'CA' in residue:
                    ca = residue['CA'].get_vector()
                    coords.append([ca[0], ca[1], ca[2]])
        return coords
