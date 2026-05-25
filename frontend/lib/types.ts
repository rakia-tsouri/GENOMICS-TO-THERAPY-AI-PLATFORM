export type Role = "researcher" | "admin";

export type JobStatus = "pending" | "running" | "completed" | "failed";

export interface User {
  id: number;
  email: string;
  full_name: string;
  role: Role;
  is_active: boolean;
  created_at: string;
}

export interface AuthToken {
  access_token: string;
  token_type: string;
}

export interface Project {
  id: number;
  name: string;
  description: string;
  cancer_type: string;
  owner_id: number;
  created_at: string;
}

export interface JobSummary {
  id: number;
  name: string;
  status: JobStatus;
  current_step: string;
  project_id: number | null;
  gene_id: string | null;
  created_at: string;
  completed_at: string | null;
}

// ----- Stage result shapes (loosely typed; backend payloads may extend) -----

export interface OrfResult {
  protein_seq: string;
  protein_length: number;
  [key: string]: unknown;
}

export interface ValidationResult {
  valid: boolean;
  gc_percent: number;
  dna_length: number;
  orfs: OrfResult[];
  warnings?: string[];
  errors?: string[];
  [key: string]: unknown;
}

export interface AnalysisResult {
  gene_id: string;
  blast: {
    protein_status: string;
    top_hit_name: string;
    identity_percent: number;
    [key: string]: unknown;
  };
  annotation: {
    uniprot_id: string;
    function: string;
    diseases: string[];
    domains: string[];
    [key: string]: unknown;
  };
  fold_check: {
    foldable: boolean;
    disordered_percent: number;
    [key: string]: unknown;
  };
  structure_3d: {
    source: string;
    pdb_file_path: string;
    plddt_mean: number;
    confidence_level: string;
    [key: string]: unknown;
  };
  [key: string]: unknown;
}

export interface DrugCandidate {
  chembl_id: string;
  smiles: string;
  binding_score: number;
  toxicity_risk: string;
  lipinski_pass: boolean;
  admet?: Record<string, unknown>;
  [key: string]: unknown;
}

export interface DrugResult {
  pocket: {
    score: number;
    volume: number;
    residues: string[];
    [key: string]: unknown;
  };
  top_candidates: DrugCandidate[];
  visualization_html_path: string;
  [key: string]: unknown;
}

export interface MutationPrediction {
  gene: string;
  mutated: boolean;
  probability: number;
  confidence: number;
  [key: string]: unknown;
}

export interface HistopathologyResult {
  model: string;
  num_patches: number;
  patches_kept: number;
  mutations: MutationPrediction[];
  gradcam_overlay_path: string;
  warnings?: string[];
  [key: string]: unknown;
}

export interface FusionGene {
  gene: string;
  genomic_signal: string;
  visual_signal: string;
  visual_probability: number;
  agreement: string;
  combined_confidence: number;
  flag_for_review: boolean;
  note: string;
  [key: string]: unknown;
}

export interface FusionResult {
  overall_agreement: number;
  cohen_kappa: number;
  flagged_genes: string[];
  genes: FusionGene[];
  [key: string]: unknown;
}

export interface Job extends JobSummary {
  error: string;
  dna_sequence: string | null;
  wsi_image_path: string | null;
  target_mutations: string[] | null;
  validation_result: ValidationResult | null;
  analysis_result: AnalysisResult | null;
  drug_result: DrugResult | null;
  histopathology_result: HistopathologyResult | null;
  fusion_result: FusionResult | null;
}

export interface Report {
  id: number;
  job_id: number;
  owner_id: number;
  title: string;
  summary: Record<string, unknown> | null;
  created_at: string;
}

// ----- Request payload shapes -----

export interface CreateProjectInput {
  name: string;
  description: string;
  cancer_type: string;
}

export interface CreateJobInput {
  name: string;
  project_id?: number | null;
  dna_sequence?: string | null;
  gene_id?: string | null;
  wsi_image_path?: string | null;
  target_mutations: string[];
}

export interface UploadResult {
  wsi_image_path: string;
  filename: string;
  size_bytes: number;
}

export interface UpdateUserInput {
  full_name?: string;
  role?: Role;
  is_active?: boolean;
}

export interface UpdateMeInput {
  full_name?: string;
  password?: string;
}
