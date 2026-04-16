import io
from Bio.Blast import NCBIXML
from Bio.PDB import PDBParser
from typing import List, Dict, Any, Optional

def parse_blast_xml(xml_data: str) -> List[Dict[str, Any]]:
    """Parse BLAST XML and return top hits."""
    result_handle = io.StringIO(xml_data)
    blast_records = NCBIXML.parse(result_handle)
    hits = []
    
    for record in blast_records:
        for alignment in record.alignments:
            hsp = alignment.hsps[0] # Top HSP
            
            # Extract organism from description (usually in brackets)
            title = alignment.title
            organism = "Unknown"
            if "[" in title and "]" in title:
                organism = title.split("[")[-1].split("]")[0]
            
            # Uniprot ID extraction (e.g. sp|P04637|P53_HUMAN)
            uniprot_id = None
            if "sp|" in title or "tr|" in title:
                parts = title.split("|")
                if len(parts) >= 2:
                    uniprot_id = parts[1]

            hits.append({
                "name": alignment.hit_def,
                "organism": organism,
                "identity_percent": (hsp.identities / hsp.align_length) * 100,
                "e_value": hsp.expect,
                "coverage_percent": (hsp.align_length / record.query_length) * 100,
                "uniprot_id": uniprot_id
            })
    return hits

def extract_plddt_from_pdb(pdb_content: str) -> Dict[str, Any]:
    """
    Extract pLDDT scores from a PDB string.
    In AlphaFold/ESMFold PDBs, the B-factor field (columns 61-66) 
    stores the pLDDT score.
    """
    parser = PDBParser(QUIET=True)
    pdb_handle = io.StringIO(pdb_content)
    structure = parser.get_structure("protein", pdb_handle)
    
    plddt_scores = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Get the B-factor of the first atom (usually CA is used, but all atoms have it)
                # We'll take the average of all atoms in the residue or just one.
                # Standard practice for AF2 is that all atoms in a residue have the same pLDDT.
                atoms = list(residue.get_atoms())
                if atoms:
                    plddt_scores.append(atoms[0].get_bfactor())
    
    if not plddt_scores:
        return {"mean": 0, "per_residue": []}
        
    mean_plddt = sum(plddt_scores) / len(plddt_scores)
    
    confidence = "low"
    if mean_plddt > 90:
        confidence = "very high"
    elif mean_plddt > 70:
        confidence = "confident"
    elif mean_plddt > 50:
        confidence = "low"
    else:
        confidence = "very low"
        
    return {
        "mean": round(mean_plddt, 2),
        "per_residue": [round(s, 2) for s in plddt_scores],
        "confidence": confidence
    }
