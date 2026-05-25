"use client";

import { Suspense, useEffect, useRef, useState } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import {
  Dna,
  Microscope,
  Upload,
  FileCheck2,
  X,
  Target,
} from "lucide-react";
import { get, post, uploadFile, ApiError } from "@/lib/api";
import type { CreateJobInput, Job, Project, UploadResult } from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Card, { CardHeader, CardBody } from "@/components/Card";
import { Field, Input, Textarea, Select } from "@/components/Input";
import Badge from "@/components/Badge";
import { ErrorBanner, InfoBanner, Spinner } from "@/components/Feedback";
import { cn, formatBytes } from "@/lib/utils";

const DEFAULT_MUTATIONS = ["TP53", "IDH1", "KRAS"];
const SUGGESTED_MUTATIONS = [
  "TP53",
  "IDH1",
  "KRAS",
  "EGFR",
  "BRAF",
  "PTEN",
  "PIK3CA",
  "ATRX",
  "NRAS",
];

type GenomicsMode = "dna" | "gene";

function NewJobForm() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const preselectedProject = searchParams.get("project_id");

  const [projects, setProjects] = useState<Project[]>([]);
  const [name, setName] = useState("");
  const [projectId, setProjectId] = useState<string>(preselectedProject ?? "");

  const [genomicsMode, setGenomicsMode] = useState<GenomicsMode>("dna");
  const [dnaSequence, setDnaSequence] = useState("");
  const [geneId, setGeneId] = useState("");

  const [upload, setUpload] = useState<UploadResult | null>(null);
  const [uploading, setUploading] = useState(false);
  const [uploadError, setUploadError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const [mutations, setMutations] = useState<string[]>(DEFAULT_MUTATIONS);
  const [customMutation, setCustomMutation] = useState("");

  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    get<Project[]>("/projects")
      .then((p) => {
        if (!cancelled) setProjects(p);
      })
      .catch(() => {
        /* non-fatal */
      });
    return () => {
      cancelled = true;
    };
  }, []);

  function toggleMutation(m: string) {
    setMutations((prev) =>
      prev.includes(m) ? prev.filter((x) => x !== m) : [...prev, m]
    );
  }

  function addCustomMutation() {
    const v = customMutation.trim().toUpperCase();
    if (v && !mutations.includes(v)) {
      setMutations((prev) => [...prev, v]);
    }
    setCustomMutation("");
  }

  async function onFileChange(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (!file) return;
    setUploadError(null);
    setUploading(true);
    try {
      const result = await uploadFile(file);
      setUpload(result);
    } catch (err) {
      setUploadError(
        err instanceof ApiError ? err.message : "Upload failed. Try again."
      );
    } finally {
      setUploading(false);
      if (fileInputRef.current) fileInputRef.current.value = "";
    }
  }

  function clearUpload() {
    setUpload(null);
    setUploadError(null);
  }

  const hasGenomics =
    genomicsMode === "dna" ? dnaSequence.trim().length > 0 : geneId.trim().length > 0;
  const hasHisto = !!upload;
  const hasAnyInput = hasGenomics || hasHisto;

  async function onSubmit(e: React.FormEvent) {
    e.preventDefault();
    setError(null);

    if (!name.trim()) {
      setError("Please give this analysis a name.");
      return;
    }
    if (!hasAnyInput) {
      setError(
        "Provide at least one input: a DNA sequence, a gene ID, or a WSI image."
      );
      return;
    }

    const payload: CreateJobInput = {
      name: name.trim(),
      project_id: projectId ? Number(projectId) : null,
      target_mutations: mutations,
    };
    if (genomicsMode === "dna" && dnaSequence.trim()) {
      payload.dna_sequence = dnaSequence.trim();
    }
    if (genomicsMode === "gene" && geneId.trim()) {
      payload.gene_id = geneId.trim();
    }
    if (upload) {
      payload.wsi_image_path = upload.wsi_image_path;
    }

    setSubmitting(true);
    try {
      const job = await post<Job>("/jobs", payload);
      router.push(`/jobs/${job.id}`);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to create analysis.");
      setSubmitting(false);
    }
  }

  return (
    <div>
      <PageHeader
        title="New analysis"
        description="Configure inputs for an end-to-end genomics-to-therapy pipeline run."
      />

      <form onSubmit={onSubmit} className="space-y-6">
        <ErrorBanner message={error} />

        <Card>
          <CardHeader title="Basics" description="Name and optional project." />
          <CardBody className="grid grid-cols-1 gap-4 sm:grid-cols-2">
            <Field label="Analysis name" htmlFor="name" required>
              <Input
                id="name"
                required
                value={name}
                onChange={(e) => setName(e.target.value)}
                placeholder="e.g. GBM-sample-001 full pipeline"
              />
            </Field>
            <Field label="Project" htmlFor="project" hint="Optional">
              <Select
                id="project"
                value={projectId}
                onChange={(e) => setProjectId(e.target.value)}
              >
                <option value="">No project</option>
                {projects.map((p) => (
                  <option key={p.id} value={p.id}>
                    {p.name}
                  </option>
                ))}
              </Select>
            </Field>
          </CardBody>
        </Card>

        <Card>
          <CardHeader
            title="Genomics input"
            description="Provide a DNA / FASTA sequence, or reference an existing gene."
            icon={<Dna className="h-4 w-4" />}
          />
          <CardBody className="space-y-4">
            <div className="inline-flex rounded-lg border border-slate-200 bg-slate-50 p-1">
              <button
                type="button"
                onClick={() => setGenomicsMode("dna")}
                className={cn(
                  "rounded-md px-3 py-1.5 text-sm font-medium transition-colors",
                  genomicsMode === "dna"
                    ? "bg-white text-brand-700 shadow-sm"
                    : "text-slate-500 hover:text-slate-700"
                )}
              >
                DNA / FASTA
              </button>
              <button
                type="button"
                onClick={() => setGenomicsMode("gene")}
                className={cn(
                  "rounded-md px-3 py-1.5 text-sm font-medium transition-colors",
                  genomicsMode === "gene"
                    ? "bg-white text-brand-700 shadow-sm"
                    : "text-slate-500 hover:text-slate-700"
                )}
              >
                Gene ID
              </button>
            </div>

            {genomicsMode === "dna" ? (
              <Field
                label="DNA sequence"
                htmlFor="dna"
                hint="Paste raw nucleotides or a FASTA record. ACGT/N accepted."
              >
                <Textarea
                  id="dna"
                  value={dnaSequence}
                  onChange={(e) => setDnaSequence(e.target.value)}
                  placeholder={">sample\nATGGCC...TAA"}
                  rows={6}
                />
              </Field>
            ) : (
              <Field
                label="Gene ID"
                htmlFor="gene"
                hint="e.g. a HGNC symbol or Ensembl/Entrez ID (TP53, ENSG00000141510)."
              >
                <Input
                  id="gene"
                  className="font-mono"
                  value={geneId}
                  onChange={(e) => setGeneId(e.target.value)}
                  placeholder="TP53"
                />
              </Field>
            )}
          </CardBody>
        </Card>

        <Card>
          <CardHeader
            title="Histopathology input (Track B)"
            description="Optional whole-slide image for mutation prediction."
            icon={<Microscope className="h-4 w-4" />}
          />
          <CardBody className="space-y-3">
            <div className="rounded-lg border border-amber-300 bg-amber-50 px-3 py-2 text-xs text-amber-800">
              <span className="font-semibold">Demo:</span> the Track&nbsp;B model is
              not trained yet — predictions are illustrative of the pipeline, not valid results.
            </div>
            <input
              ref={fileInputRef}
              type="file"
              accept=".svs,.tiff,.tif,.png,.jpg,.jpeg,.ndpi,image/*"
              className="hidden"
              onChange={onFileChange}
            />

            {!upload ? (
              <button
                type="button"
                onClick={() => fileInputRef.current?.click()}
                disabled={uploading}
                className="flex w-full flex-col items-center justify-center gap-2 rounded-lg border-2 border-dashed border-slate-300 bg-slate-50 px-6 py-10 text-center transition-colors hover:border-brand-300 hover:bg-brand-50/40 disabled:opacity-70"
              >
                {uploading ? (
                  <>
                    <Spinner />
                    <span className="text-sm text-slate-500">Uploading…</span>
                  </>
                ) : (
                  <>
                    <span className="flex h-10 w-10 items-center justify-center rounded-full bg-brand-100 text-brand-600">
                      <Upload className="h-5 w-5" />
                    </span>
                    <span className="text-sm font-medium text-slate-700">
                      Click to upload a WSI
                    </span>
                    <span className="text-xs text-slate-400">
                      .svs, .tiff, .ndpi or standard image formats
                    </span>
                  </>
                )}
              </button>
            ) : (
              <div className="flex items-center justify-between rounded-lg border border-green-200 bg-green-50 px-4 py-3">
                <div className="flex items-center gap-3">
                  <span className="flex h-9 w-9 items-center justify-center rounded-lg bg-green-100 text-green-700">
                    <FileCheck2 className="h-5 w-5" />
                  </span>
                  <div>
                    <p className="text-sm font-medium text-slate-900">
                      {upload.filename}
                    </p>
                    <p className="text-xs text-slate-500">
                      {formatBytes(upload.size_bytes)} ·{" "}
                      <span className="font-mono">{upload.wsi_image_path}</span>
                    </p>
                  </div>
                </div>
                <button
                  type="button"
                  onClick={clearUpload}
                  className="rounded-md p-1.5 text-slate-400 hover:bg-white hover:text-slate-600"
                  aria-label="Remove file"
                >
                  <X className="h-4 w-4" />
                </button>
              </div>
            )}

            <ErrorBanner message={uploadError} />
          </CardBody>
        </Card>

        <Card>
          <CardHeader
            title="Target mutations"
            description="Genes to evaluate across genomic and visual modalities."
            icon={<Target className="h-4 w-4" />}
          />
          <CardBody className="space-y-4">
            <div className="flex flex-wrap gap-2">
              {SUGGESTED_MUTATIONS.map((m) => {
                const active = mutations.includes(m);
                return (
                  <button
                    key={m}
                    type="button"
                    onClick={() => toggleMutation(m)}
                    className={cn(
                      "rounded-full border px-3 py-1 text-sm font-medium transition-colors",
                      active
                        ? "border-brand-300 bg-brand-50 text-brand-700"
                        : "border-slate-200 bg-white text-slate-600 hover:bg-slate-50"
                    )}
                  >
                    {m}
                  </button>
                );
              })}
            </div>

            {mutations.filter((m) => !SUGGESTED_MUTATIONS.includes(m)).length >
              0 && (
              <div className="flex flex-wrap gap-2">
                {mutations
                  .filter((m) => !SUGGESTED_MUTATIONS.includes(m))
                  .map((m) => (
                    <Badge key={m} tone="brand">
                      {m}
                      <button
                        type="button"
                        onClick={() => toggleMutation(m)}
                        className="ml-1 text-brand-500 hover:text-brand-700"
                        aria-label={`Remove ${m}`}
                      >
                        <X className="h-3 w-3" />
                      </button>
                    </Badge>
                  ))}
              </div>
            )}

            <div className="flex gap-2">
              <Input
                value={customMutation}
                onChange={(e) => setCustomMutation(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") {
                    e.preventDefault();
                    addCustomMutation();
                  }
                }}
                placeholder="Add a custom gene symbol…"
                className="max-w-xs font-mono"
              />
              <Button
                type="button"
                variant="outline"
                onClick={addCustomMutation}
                disabled={!customMutation.trim()}
              >
                Add
              </Button>
            </div>
            <p className="text-xs text-slate-500">
              {mutations.length} gene{mutations.length === 1 ? "" : "s"} selected
            </p>
          </CardBody>
        </Card>

        {!hasAnyInput && (
          <InfoBanner>
            Provide at least one input — a DNA sequence, a gene ID, or a
            whole-slide image — to launch the pipeline.
          </InfoBanner>
        )}

        <div className="flex items-center justify-end gap-3">
          <Button
            type="button"
            variant="outline"
            onClick={() => router.back()}
          >
            Cancel
          </Button>
          <Button type="submit" loading={submitting} disabled={!hasAnyInput}>
            Launch analysis
          </Button>
        </div>
      </form>
    </div>
  );
}

export default function NewJobPage() {
  return (
    <Suspense fallback={<div className="py-10"><Spinner /></div>}>
      <NewJobForm />
    </Suspense>
  );
}
