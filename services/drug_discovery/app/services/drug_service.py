import httpx
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import List, Dict, Optional
import os

class DrugService:
    CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
    ZINC_BASE = "https://zinc15.docking.org/substances/subsets"

    @staticmethod
    async def get_candidates(uniprot_id: Optional[str] = None) -> List[Dict]:
        """Step B: Fetch drug candidates from ChEMBL and ZINC15."""
        candidates = []
        
        # 1. ChEMBL: Known drugs for this protein
        if uniprot_id:
            chembl_drugs = await DrugService._fetch_chembl_drugs(uniprot_id)
            candidates.extend(chembl_drugs)
        
        # 2. ZINC15: Drug-like molecules (FDA approved subset)
        zinc_drugs = await DrugService._fetch_zinc_drugs()
        candidates.extend(zinc_drugs)
        
        # 3. Filter by Lipinski's Rule of 5
        filtered = [c for c in candidates if DrugService._lipinski_filter(c["smiles"])]
        
        return filtered

    @staticmethod
    async def _fetch_chembl_drugs(uniprot_id: str) -> List[Dict]:
        async with httpx.AsyncClient() as client:
            # UniProt ID -> ChEMBL ID
            target_url = f"{DrugService.CHEMBL_BASE}/target.json?target_components__accession={uniprot_id}"
            resp = await client.get(target_url)
            targets = resp.json().get("targets", [])
            if not targets:
                return []
            
            target_id = targets[0]["target_chembl_id"]
            
            # Activities for this target
            act_url = f"{DrugService.CHEMBL_BASE}/activity.json?target_chembl_id={target_id}&limit=50"
            resp = await client.get(act_url)
            activities = resp.json().get("activities", [])
            
            drugs = []
            for act in activities:
                pchembl = act.get("pchembl_value")
                if pchembl and 4 <= float(pchembl) <= 12:
                    drugs.append({
                        "smiles": act["canonical_smiles"],
                        "chembl_id": act["molecule_chembl_id"],
                        "pchembl_value": float(pchembl)
                    })
            return drugs

    @staticmethod
    async def _fetch_zinc_drugs() -> List[Dict]:
        async with httpx.AsyncClient() as client:
            url = f"{DrugService.ZINC_BASE}/fda.json?count=50"
            try:
                resp = await client.get(url)
                data = resp.json()
                return [{"smiles": d["smiles"], "chembl_id": None, "pchembl_value": None} for d in data]
            except:
                return []

    @staticmethod
    def _lipinski_filter(smiles: str) -> bool:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return False
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        return violations <= 1
