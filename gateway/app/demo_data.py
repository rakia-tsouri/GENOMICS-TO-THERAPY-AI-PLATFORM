"""Hand-curated demo data: realistic projects, completed analyses, and reports.

Inserted on startup by ``seed.py`` so the platform looks populated for product
demos and first-run UX. Idempotent (matched by name).
"""

# ---------------------------------------------------------------------------
# Projects (owned by the demo researcher)
# ---------------------------------------------------------------------------
PROJECTS = [
    {
        "name": "TP53 — Lung Adenocarcinoma",
        "description": (
            "Reference TP53 screen on a LUAD cohort. Targets the canonical "
            "tumor-suppressor mutation hotspots (R175, R248, R273, R282)."
        ),
        "cancer_type": "LUAD",
    },
    {
        "name": "BRCA1/2 — Breast Cancer Screen",
        "description": (
            "BRCA1/BRCA2 germline + somatic screen across a hereditary breast "
            "cancer cohort. Pairs sequence analysis with H&E histology."
        ),
        "cancer_type": "BRCA",
    },
    {
        "name": "EGFR — NSCLC Targeted Therapy",
        "description": (
            "EGFR-driven non-small-cell lung cancer cohort. Prioritises "
            "ATP-competitive and covalent EGFR inhibitors."
        ),
        "cancer_type": "NSCLC",
    },
    {
        "name": "IDH1 — Glioma Subtyping",
        "description": (
            "IDH1 R132H/R132C subtyping for diffuse glioma. Cross-validates "
            "sequence variants against tissue morphology."
        ),
        "cancer_type": "GBM",
    },
]

# ---------------------------------------------------------------------------
# Realistic protein sequence stubs (truncated; UI shows length only)
# ---------------------------------------------------------------------------
_TP53 = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTP"
    "AAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYK"
    "QSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGM"
    "NRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYF"
    "TLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)
_BRCA1 = (
    "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQL"
    "VEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRT"
    "LRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDL"
    "NTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETC"
)
_EGFR = (
    "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTI"
    "QEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQW"
    "RDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRES"
)
_IDH1 = (
    "MSKKISGGSVVEMQGDEMTRIIWELIKEKLIFPYVELDLHSYDLGIENRDATNDQVTKDAAEAIKKHNVGVKCATITPDEKR"
    "VEEFKLKQMWKSPNGTIRNILGGTVFREAIICKNIPRLVSGWVKPIIIGRHAYGDQYRATDFVVPGPGKVEITYTPSDGTQK"
    "VTYLVHNFEEGGGVAMGMYNQDKSIEDFAHSSFQMALSKGWPLYLSTKNTILKKYDGRFKDIFQEIYDKQYKSQFEAQKIWY"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _validation(*, gene_id, dna_length, gc, protein):
    return {
        "valid": True,
        "gene_id": gene_id,
        "dna_length": dna_length,
        "gc_percent": gc,
        "orfs": [{
            "start": 0,
            "end": dna_length,
            "length": dna_length,
            "dna_seq": f"ATG{'NNN' * 6}...{'NNN' * 6}TAA",
            "protein_seq": protein,
            "protein_length": len(protein),
        }],
        "warnings": [],
        "errors": [],
    }


def _analysis(*, gene_id, protein, blast, annotation, plddt, conf, source="AlphaFold DB", pdb_path):
    return {
        "gene_id": gene_id,
        "protein_seq": protein,
        "protein_length": len(protein),
        "blast": blast,
        "annotation": annotation,
        "fold_check": {
            "foldindex_score": 0.62,
            "foldable": True,
            "disordered_percent": 18.4,
            "disordered_regions": [[1, 35]],
            "method": "api",
            "warning": None,
        },
        "structure_3d": {
            "source": source,
            "status": "success",
            "uniprot_id": annotation.get("uniprot_id"),
            "pdb_file_path": pdb_path,
            "plddt_mean": plddt,
            "plddt_per_residue": [],
            "confidence_level": conf,
            "reason": None,
            "fallback": False,
        },
        "processing_time_sec": 14.2,
    }


def _candidate(chembl_id, smiles, binding, pchembl, mw, logp, tox, lipinski=True, bioav=0.78):
    return {
        "chembl_id": chembl_id,
        "smiles": smiles,
        "binding_score": binding,
        "pchembl_value": pchembl,
        "lipinski_pass": lipinski,
        "admet": {
            "mw": mw, "logp": logp, "hbd": 2, "hba": 6,
            "tpsa": 75.4, "bioavailability": bioav,
        },
        "toxicity_risk": tox,
    }


def _drug(*, gene_id, pocket_residues, score, volume, candidates, viz):
    return {
        "gene_id": gene_id,
        "pocket": {
            "residues": pocket_residues,
            "center": [12.4, -3.7, 8.1],
            "volume": volume,
            "score": score,
        },
        "top_candidates": candidates,
        "visualization_html_path": viz,
    }


def _histo(*, num, kept, mutations, model_tag, gradcam):
    return {
        "status": "success",
        "num_patches": num,
        "patches_kept": kept,
        "mutations": mutations,
        "gradcam_overlay_path": gradcam,
        "model": model_tag,
        "processing_time_sec": 21.4,
        "warnings": [],
        "reason": None,
    }


def _fusion(*, genes, overall, kappa, flagged):
    return {
        "genes": genes,
        "overall_agreement": overall,
        "cohen_kappa": kappa,
        "flagged_genes": flagged,
    }


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------
JOBS = [
    # ------------------------- 1 — TP53 / LUAD (concordant) -------------------------
    {
        "project": "TP53 — Lung Adenocarcinoma",
        "days_ago": 6,
        "name": "TP53 reference analysis · LUAD",
        "gene_id": "NM_000546.6",
        "wsi_image_path": "/data/structures/uploads/luad_tp53_demo.png",
        "target_mutations": ["TP53", "KRAS", "EGFR"],
        "validation_result": _validation(gene_id="NM_000546.6", dna_length=2591, gc=51.2, protein=_TP53),
        "analysis_result": _analysis(
            gene_id="NM_000546.6", protein=_TP53,
            blast={
                "top_hit_name": "Cellular tumor antigen p53",
                "top_hit_organism": "Homo sapiens",
                "identity_percent": 100.0,
                "e_value": 0.0,
                "coverage_percent": 100.0,
                "uniprot_id": "P04637",
                "protein_status": "known",
            },
            annotation={
                "uniprot_id": "P04637",
                "function": "Tumor suppressor; induces cell-cycle arrest or apoptosis in response to DNA damage.",
                "domains": ["Transactivation", "DNA-binding domain", "Tetramerization", "Regulatory"],
                "active_sites": [{"position": 175, "description": "DNA-binding interface"}],
                "binding_sites": [{"positions": [248, 273, 282], "ligand": "Zn2+"}],
                "diseases": ["Li-Fraumeni syndrome 1", "Lung adenocarcinoma", "Colorectal cancer",
                             "Hepatocellular carcinoma", "Breast cancer"],
            },
            plddt=91.2, conf="very high", pdb_path="/demo/structures/1TUP.pdb",
        ),
        "drug_result": _drug(
            gene_id="NM_000546.6", pocket_residues=[175, 176, 245, 248, 273, 282], score=0.94, volume=842.3,
            candidates=[
                _candidate("CHEMBL428647",
                           "CC(C)c1nc(c(n1[C@@H]1Cc2ccc(Cl)cc2-c2cc(Cl)ccc12)C(=O)N1CCC(O)CC1)C",
                           7.42, 6.81, 581.5, 4.6, "low"),
                _candidate("CHEMBL3989958",
                           "CC1(C(=O)N(C2CC2)C3(CC3)C(=O)N1C4=NC(=NC=C4Cl)N5CCC(CC5)(F)F)C6=CC=C(C=C6)C#N",
                           6.95, 6.12, 582.6, 3.2, "low"),
                _candidate("CHEMBL1336",
                           "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",
                           6.71, 5.94, 464.8, 4.3, "moderate"),
                _candidate("CHEMBL941",
                           "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
                           6.28, 5.71, 493.6, 3.0, "low"),
                _candidate("CHEMBL553",
                           "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
                           5.89, 5.42, 393.4, 3.4, "low"),
            ],
            viz="/data/structures/visualizations/viz_tp53_luad.html",
        ),
        "histopathology_result": _histo(
            num=184, kept=156,
            mutations=[
                {"gene": "TP53", "mutated": True,  "probability": 0.92, "confidence": "high"},
                {"gene": "KRAS", "mutated": False, "probability": 0.21, "confidence": "high"},
                {"gene": "EGFR", "mutated": False, "probability": 0.18, "confidence": "high"},
            ],
            model_tag="ResNet50 + Attention-MIL (TCGA-LUAD)",
            gradcam="/demo/gradcam/gradcam_tp53_luad.png",
        ),
        "fusion_result": _fusion(
            genes=[
                {"gene": "TP53", "genomic_signal": True, "visual_signal": True, "visual_probability": 0.92,
                 "agreement": "concordant", "combined_confidence": 0.95, "flag_for_review": False,
                 "note": "Both modalities agree."},
                {"gene": "KRAS", "genomic_signal": None, "visual_signal": False, "visual_probability": 0.21,
                 "agreement": "single_modality", "combined_confidence": 0.85, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
                {"gene": "EGFR", "genomic_signal": None, "visual_signal": False, "visual_probability": 0.18,
                 "agreement": "single_modality", "combined_confidence": 0.88, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
            ],
            overall="concordant", kappa=1.0, flagged=[],
        ),
    },

    # ------------------------- 2 — TP53 / LUAD (patient sample, discordant) -------------------------
    {
        "project": "TP53 — Lung Adenocarcinoma",
        "days_ago": 4,
        "name": "TP53 patient sample · LUAD-A3F1",
        "gene_id": "NM_000546.6",
        "wsi_image_path": "/data/structures/uploads/luad_a3f1.png",
        "target_mutations": ["TP53", "KRAS"],
        "validation_result": _validation(gene_id="NM_000546.6", dna_length=2591, gc=51.4, protein=_TP53),
        "analysis_result": _analysis(
            gene_id="NM_000546.6", protein=_TP53,
            blast={
                "top_hit_name": "Cellular tumor antigen p53",
                "top_hit_organism": "Homo sapiens",
                "identity_percent": 99.5,
                "e_value": 0.0,
                "coverage_percent": 100.0,
                "uniprot_id": "P04637",
                "protein_status": "known",
            },
            annotation={
                "uniprot_id": "P04637",
                "function": "Tumor suppressor; induces cell-cycle arrest or apoptosis in response to DNA damage.",
                "domains": ["Transactivation", "DNA-binding domain", "Tetramerization"],
                "active_sites": [],
                "binding_sites": [{"positions": [248, 273], "ligand": "Zn2+"}],
                "diseases": ["Li-Fraumeni syndrome 1", "Lung adenocarcinoma"],
            },
            plddt=88.7, conf="confident", pdb_path="/demo/structures/1TUP.pdb",
        ),
        "drug_result": _drug(
            gene_id="NM_000546.6", pocket_residues=[175, 248, 273, 282], score=0.88, volume=794.1,
            candidates=[
                _candidate("CHEMBL428647",
                           "CC(C)c1nc(c(n1[C@@H]1Cc2ccc(Cl)cc2-c2cc(Cl)ccc12)C(=O)N1CCC(O)CC1)C",
                           7.18, 6.62, 581.5, 4.6, "low"),
                _candidate("CHEMBL3989958",
                           "CC1(C(=O)N(C2CC2)C3(CC3)C(=O)N1C4=NC(=NC=C4Cl)N5CCC(CC5)(F)F)C6=CC=C(C=C6)C#N",
                           6.72, 5.98, 582.6, 3.2, "low"),
                _candidate("CHEMBL553",
                           "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
                           5.61, 5.20, 393.4, 3.4, "low"),
            ],
            viz="/data/structures/visualizations/viz_tp53_a3f1.html",
        ),
        "histopathology_result": _histo(
            num=212, kept=178,
            mutations=[
                {"gene": "TP53", "mutated": False, "probability": 0.27, "confidence": "medium"},
                {"gene": "KRAS", "mutated": True,  "probability": 0.81, "confidence": "high"},
            ],
            model_tag="ResNet50 + Attention-MIL (TCGA-LUAD)",
            gradcam="/demo/gradcam/gradcam_a3f1.png",
        ),
        "fusion_result": _fusion(
            genes=[
                {"gene": "TP53", "genomic_signal": True, "visual_signal": False, "visual_probability": 0.27,
                 "agreement": "discordant", "combined_confidence": 0.35, "flag_for_review": True,
                 "note": "Genomic and visual evidence disagree — manual review recommended."},
                {"gene": "KRAS", "genomic_signal": None, "visual_signal": True, "visual_probability": 0.81,
                 "agreement": "single_modality", "combined_confidence": 0.78, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
            ],
            overall="discordant", kappa=0.0, flagged=["TP53"],
        ),
    },

    # ------------------------- 3 — BRCA1 / breast cancer -------------------------
    {
        "project": "BRCA1/2 — Breast Cancer Screen",
        "days_ago": 3,
        "name": "BRCA1 hereditary screen · BC-072",
        "gene_id": "NM_007294.4",
        "wsi_image_path": "/data/structures/uploads/brca_072.png",
        "target_mutations": ["BRCA1", "BRCA2", "TP53"],
        "validation_result": _validation(gene_id="NM_007294.4", dna_length=7224, gc=41.2, protein=_BRCA1),
        "analysis_result": _analysis(
            gene_id="NM_007294.4", protein=_BRCA1,
            blast={
                "top_hit_name": "Breast cancer type 1 susceptibility protein",
                "top_hit_organism": "Homo sapiens",
                "identity_percent": 99.9,
                "e_value": 0.0,
                "coverage_percent": 100.0,
                "uniprot_id": "P38398",
                "protein_status": "known",
            },
            annotation={
                "uniprot_id": "P38398",
                "function": "E3 ubiquitin-protein ligase; central in homologous-recombination DNA repair.",
                "domains": ["RING-type zinc finger", "BRCT 1", "BRCT 2"],
                "active_sites": [{"position": 24, "description": "RING ligase catalytic"}],
                "binding_sites": [{"positions": [1656, 1700], "ligand": "phospho-peptide"}],
                "diseases": ["Hereditary breast and ovarian cancer", "Fanconi anemia complementation group S",
                             "Pancreatic cancer"],
            },
            plddt=87.9, conf="confident", pdb_path="/demo/structures/1JM7.pdb",
        ),
        "drug_result": _drug(
            gene_id="NM_007294.4", pocket_residues=[1656, 1657, 1700, 1701, 1739, 1740], score=0.86, volume=712.6,
            candidates=[
                _candidate("CHEMBL521686",
                           "C1CC1C(=O)N2CCN(CC2)Cc3ccc(F)c(C4=NNC(=O)c5ccccc54)c3",
                           7.05, 6.74, 434.5, 1.9, "low"),
                _candidate("CHEMBL2105760",
                           "Fc1ccc(F)c(C2=NNc3cc(NC(=O)C4CC4)ccc23)c1",
                           6.62, 6.10, 365.4, 2.8, "low"),
                _candidate("CHEMBL3989523",
                           "O=C1NN=C(c2ccccc12)c1ccc(N2CCC(N3CCN(C)CC3)CC2)cc1F",
                           6.34, 5.81, 421.5, 2.6, "low"),
                _candidate("CHEMBL83",
                           "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
                           5.78, 4.92, 371.5, 6.7, "moderate", lipinski=False),
            ],
            viz="/data/structures/visualizations/viz_brca1_072.html",
        ),
        "histopathology_result": _histo(
            num=168, kept=142,
            mutations=[
                {"gene": "BRCA1", "mutated": True,  "probability": 0.87, "confidence": "high"},
                {"gene": "BRCA2", "mutated": False, "probability": 0.38, "confidence": "medium"},
                {"gene": "TP53",  "mutated": False, "probability": 0.29, "confidence": "medium"},
            ],
            model_tag="ResNet50 + Attention-MIL (TCGA-BRCA)",
            gradcam="/demo/gradcam/gradcam_brca_072.png",
        ),
        "fusion_result": _fusion(
            genes=[
                {"gene": "BRCA1", "genomic_signal": True, "visual_signal": True, "visual_probability": 0.87,
                 "agreement": "concordant", "combined_confidence": 0.93, "flag_for_review": False,
                 "note": "Both modalities agree."},
                {"gene": "BRCA2", "genomic_signal": None, "visual_signal": False, "visual_probability": 0.38,
                 "agreement": "single_modality", "combined_confidence": 0.62, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
                {"gene": "TP53", "genomic_signal": None, "visual_signal": False, "visual_probability": 0.29,
                 "agreement": "single_modality", "combined_confidence": 0.70, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
            ],
            overall="concordant", kappa=1.0, flagged=[],
        ),
    },

    # ------------------------- 4 — EGFR / NSCLC -------------------------
    {
        "project": "EGFR — NSCLC Targeted Therapy",
        "days_ago": 2,
        "name": "EGFR L858R · NSCLC-114",
        "gene_id": "NM_005228.5",
        "wsi_image_path": "/data/structures/uploads/nsclc_114.png",
        "target_mutations": ["EGFR", "KRAS"],
        "validation_result": _validation(gene_id="NM_005228.5", dna_length=5616, gc=58.7, protein=_EGFR),
        "analysis_result": _analysis(
            gene_id="NM_005228.5", protein=_EGFR,
            blast={
                "top_hit_name": "Epidermal growth factor receptor",
                "top_hit_organism": "Homo sapiens",
                "identity_percent": 99.7,
                "e_value": 0.0,
                "coverage_percent": 100.0,
                "uniprot_id": "P00533",
                "protein_status": "known",
            },
            annotation={
                "uniprot_id": "P00533",
                "function": "Receptor tyrosine kinase binding EGF; drives proliferation and survival.",
                "domains": ["Receptor L-domain 1", "Furin-like cysteine-rich", "Tyrosine kinase domain"],
                "active_sites": [{"position": 745, "description": "ATP binding"}],
                "binding_sites": [{"positions": [718, 745, 858], "ligand": "ATP"}],
                "diseases": ["Non-small cell lung cancer", "Adenocarcinoma of lung",
                             "Glioblastoma", "Inflammatory skin and bowel disease"],
            },
            plddt=89.4, conf="confident", pdb_path="/demo/structures/1M17.pdb",
        ),
        "drug_result": _drug(
            gene_id="NM_005228.5", pocket_residues=[718, 719, 745, 790, 858, 859], score=0.92, volume=915.8,
            candidates=[
                _candidate("CHEMBL3353410",
                           "CN(C)CCN(C)c1c(NC(=O)C=C)cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)c1",
                           7.84, 7.42, 499.6, 3.6, "low"),
                _candidate("CHEMBL939",
                           "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
                           7.31, 6.95, 446.9, 4.1, "low"),
                _candidate("CHEMBL553",
                           "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
                           6.92, 6.41, 393.4, 3.4, "low"),
                _candidate("CHEMBL428690",
                           "Cc1ccc(C(=O)Nc2cc3c(Nc4ccc(F)c(Cl)c4)ncnc3cc2OCCCN2CCOCC2)cc1",
                           6.58, 6.04, 581.1, 4.7, "low"),
            ],
            viz="/data/structures/visualizations/viz_egfr_114.html",
        ),
        "histopathology_result": _histo(
            num=196, kept=164,
            mutations=[
                {"gene": "EGFR", "mutated": True,  "probability": 0.91, "confidence": "high"},
                {"gene": "KRAS", "mutated": False, "probability": 0.19, "confidence": "high"},
            ],
            model_tag="ResNet50 + Attention-MIL (TCGA-LUAD)",
            gradcam="/demo/gradcam/gradcam_nsclc_114.png",
        ),
        "fusion_result": _fusion(
            genes=[
                {"gene": "EGFR", "genomic_signal": True, "visual_signal": True, "visual_probability": 0.91,
                 "agreement": "concordant", "combined_confidence": 0.95, "flag_for_review": False,
                 "note": "Both modalities agree."},
                {"gene": "KRAS", "genomic_signal": None, "visual_signal": False, "visual_probability": 0.19,
                 "agreement": "single_modality", "combined_confidence": 0.87, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
            ],
            overall="concordant", kappa=1.0, flagged=[],
        ),
    },

    # ------------------------- 5 — IDH1 / Glioma -------------------------
    {
        "project": "IDH1 — Glioma Subtyping",
        "days_ago": 1,
        "name": "IDH1 R132H · glioma case GB-007",
        "gene_id": "NM_005896.4",
        "wsi_image_path": "/data/structures/uploads/glioma_007.png",
        "target_mutations": ["IDH1", "TP53"],
        "validation_result": _validation(gene_id="NM_005896.4", dna_length=2655, gc=46.8, protein=_IDH1),
        "analysis_result": _analysis(
            gene_id="NM_005896.4", protein=_IDH1,
            blast={
                "top_hit_name": "Isocitrate dehydrogenase [NADP] cytoplasmic",
                "top_hit_organism": "Homo sapiens",
                "identity_percent": 99.8,
                "e_value": 0.0,
                "coverage_percent": 100.0,
                "uniprot_id": "O75874",
                "protein_status": "known",
            },
            annotation={
                "uniprot_id": "O75874",
                "function": "Catalyses oxidative decarboxylation of isocitrate; R132 mutations confer neomorphic 2-HG production.",
                "domains": ["Isocitrate/isopropylmalate dehydrogenase"],
                "active_sites": [{"position": 132, "description": "Isocitrate binding (R132 hotspot)"}],
                "binding_sites": [{"positions": [100, 109, 132, 212], "ligand": "isocitrate / NADP+"}],
                "diseases": ["Glioma", "Acute myeloid leukemia", "Chondrosarcoma", "Cholangiocarcinoma"],
            },
            plddt=92.6, conf="very high", pdb_path="/demo/structures/1T0L.pdb",
        ),
        "drug_result": _drug(
            gene_id="NM_005896.4", pocket_residues=[100, 109, 132, 133, 212, 275], score=0.95, volume=1024.5,
            candidates=[
                _candidate("CHEMBL3989958",
                           "CC1(C(=O)N(C2CC2)C3(CC3)C(=O)N1C4=NC(=NC=C4Cl)N5CCC(CC5)(F)F)C6=CC=C(C=C6)C#N",
                           8.12, 7.65, 582.6, 3.2, "low"),
                _candidate("CHEMBL4297468",
                           "CC1(C2=NC=CC(=N2)C3=NC4=C(C=N3)N(C=C4Cl)C5CCNCC5)C6=CC=C(C=C6F)Cl",
                           7.61, 7.12, 502.4, 3.7, "low"),
                _candidate("CHEMBL2105845",
                           "Cc1nc(C(=O)Nc2cccc(C(F)(F)F)c2)c(c1)c1ccc(F)cc1Cl",
                           6.94, 6.32, 437.8, 4.2, "low"),
                _candidate("CHEMBL553",
                           "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
                           5.73, 5.10, 393.4, 3.4, "low"),
            ],
            viz="/data/structures/visualizations/viz_idh1_007.html",
        ),
        "histopathology_result": _histo(
            num=152, kept=131,
            mutations=[
                {"gene": "IDH1", "mutated": True,  "probability": 0.89, "confidence": "high"},
                {"gene": "TP53", "mutated": False, "probability": 0.31, "confidence": "medium"},
            ],
            model_tag="ResNet50 + Attention-MIL (TCGA-LGG)",
            gradcam="/demo/gradcam/gradcam_glioma_007.png",
        ),
        "fusion_result": _fusion(
            genes=[
                {"gene": "IDH1", "genomic_signal": True, "visual_signal": True, "visual_probability": 0.89,
                 "agreement": "concordant", "combined_confidence": 0.96, "flag_for_review": False,
                 "note": "Both modalities agree."},
                {"gene": "TP53", "genomic_signal": None, "visual_signal": False, "visual_probability": 0.31,
                 "agreement": "single_modality", "combined_confidence": 0.69, "flag_for_review": False,
                 "note": "Visual (Track B) evidence only."},
            ],
            overall="concordant", kappa=1.0, flagged=[],
        ),
    },
]

# Jobs whose names should also have a generated report.
REPORT_FOR_JOBS = [
    "TP53 reference analysis · LUAD",
    "BRCA1 hereditary screen · BC-072",
    "EGFR L858R · NSCLC-114",
]
