"use client";

import { useCallback, useEffect, useRef, useState } from "react";
import Link from "next/link";
import { useParams, useRouter } from "next/navigation";
import {
  ArrowLeft,
  Dna,
  Microscope,
  Pill,
  Boxes,
  GitCompareArrows,
  FileText,
  Trash2,
  CircleDot,
  Layers,
  ExternalLink,
} from "lucide-react";
import { get, post, del, ApiError } from "@/lib/api";
import type {
  Job,
  Report,
} from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Card, { CardHeader, CardBody } from "@/components/Card";
import Badge, { toxicityTone, agreementTone } from "@/components/Badge";
import StatusBadge from "@/components/StatusBadge";
import Stepper from "@/components/Stepper";
import ConfidenceBar from "@/components/ConfidenceBar";
import Modal from "@/components/Modal";
import { Table, THead, TBody, TR, TH, TD, EmptyRow } from "@/components/Table";
import {
  ErrorBanner,
  InfoBanner,
  PageSpinner,
  Spinner,
} from "@/components/Feedback";
import {
  basename,
  cn,
  formatDate,
  formatPercent,
  truncateMiddle,
} from "@/lib/utils";

const POLL_INTERVAL = 3000;

function KV({ label, value }: { label: string; value: React.ReactNode }) {
  return (
    <div className="flex items-baseline justify-between gap-4 py-1.5">
      <dt className="text-xs font-medium uppercase tracking-wide text-slate-500">
        {label}
      </dt>
      <dd className="text-right text-sm font-medium text-slate-900">{value}</dd>
    </div>
  );
}

export default function JobDetailPage() {
  const params = useParams<{ id: string }>();
  const router = useRouter();
  const jobId = Number(params.id);

  const [job, setJob] = useState<Job | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const [generating, setGenerating] = useState(false);
  const [reportError, setReportError] = useState<string | null>(null);
  const [confirmDelete, setConfirmDelete] = useState(false);
  const [deleting, setDeleting] = useState(false);

  const timerRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const fetchJob = useCallback(async () => {
    try {
      const data = await get<Job>(`/jobs/${jobId}`);
      setJob(data);
      setError(null);
      return data;
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to load analysis.");
      return null;
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  useEffect(() => {
    if (Number.isNaN(jobId)) return;
    let active = true;

    async function tick() {
      const data = await fetchJob();
      if (!active) return;
      if (data && (data.status === "pending" || data.status === "running")) {
        timerRef.current = setTimeout(tick, POLL_INTERVAL);
      }
    }

    void tick();

    return () => {
      active = false;
      if (timerRef.current) clearTimeout(timerRef.current);
    };
  }, [jobId, fetchJob]);

  async function onGenerateReport() {
    if (!job) return;
    setGenerating(true);
    setReportError(null);
    try {
      const report = await post<Report>("/reports", { job_id: job.id });
      router.push(`/reports/${report.id}`);
    } catch (err) {
      setReportError(
        err instanceof ApiError ? err.message : "Failed to generate report."
      );
      setGenerating(false);
    }
  }

  async function onDelete() {
    if (!job) return;
    setDeleting(true);
    try {
      await del(`/jobs/${job.id}`);
      router.replace("/jobs");
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to delete analysis.");
      setDeleting(false);
    }
  }

  if (loading) return <PageSpinner label="Loading analysis…" />;

  if (!job) {
    return (
      <div>
        <ErrorBanner message={error || "Analysis not found."} />
        <Link
          href="/jobs"
          className="mt-4 inline-flex items-center gap-1 text-sm text-brand-600"
        >
          <ArrowLeft className="h-4 w-4" /> Back to analyses
        </Link>
      </div>
    );
  }

  const isRunning = job.status === "pending" || job.status === "running";
  const isCompleted = job.status === "completed";
  const isFailed = job.status === "failed";

  return (
    <div>
      <Link
        href="/jobs"
        className="mb-4 inline-flex items-center gap-1 text-sm text-slate-500 hover:text-slate-700"
      >
        <ArrowLeft className="h-4 w-4" /> Analyses
      </Link>

      <PageHeader
        title={
          <span className="flex flex-wrap items-center gap-3">
            {job.name}
            <StatusBadge status={job.status} />
          </span>
        }
        description={
          <span className="flex flex-wrap items-center gap-x-4 gap-y-1">
            <span>Created {formatDate(job.created_at)}</span>
            {job.gene_id && (
              <span className="font-mono text-xs">Gene: {job.gene_id}</span>
            )}
            {job.completed_at && (
              <span>Completed {formatDate(job.completed_at)}</span>
            )}
          </span>
        }
        action={
          <div className="flex gap-2">
            {isCompleted && (
              <Button onClick={onGenerateReport} loading={generating}>
                <FileText className="h-4 w-4" /> Generate report
              </Button>
            )}
            <Button variant="danger" onClick={() => setConfirmDelete(true)}>
              <Trash2 className="h-4 w-4" /> Delete
            </Button>
          </div>
        }
      />

      <ErrorBanner message={reportError} className="mb-4" />

      {/* Pipeline progress */}
      <Card className="mb-6">
        <CardBody>
          <div className="mb-4 flex items-center justify-between">
            <h3 className="text-sm font-semibold text-slate-900">
              Pipeline progress
            </h3>
            {isRunning && (
              <span className="flex items-center gap-2 text-xs text-blue-600">
                <Spinner className="h-4 w-4" />
                Live · current step:{" "}
                <span className="font-mono">{job.current_step || "queued"}</span>
              </span>
            )}
          </div>
          <Stepper currentStep={job.current_step} status={job.status} />
        </CardBody>
      </Card>

      {isFailed && (
        <ErrorBanner
          className="mb-6"
          message={`Analysis failed${job.error ? `: ${job.error}` : "."}`}
        />
      )}

      {isRunning && (
        <InfoBanner className="mb-6">
          This analysis is still running. Results will appear automatically as
          each stage completes (refreshing every {POLL_INTERVAL / 1000}s).
        </InfoBanner>
      )}

      {/* Result cards */}
      <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
        {/* Genomics */}
        {job.validation_result && (
          <Card>
            <CardHeader
              title="Genomics & validation"
              icon={<Dna className="h-4 w-4" />}
            />
            <CardBody>
              <dl className="divide-y divide-slate-100">
                <KV
                  label="Sequence valid"
                  value={
                    <Badge tone={job.validation_result.valid ? "green" : "red"}>
                      {job.validation_result.valid ? "Valid" : "Invalid"}
                    </Badge>
                  }
                />
                <KV
                  label="GC content"
                  value={formatPercent(job.validation_result.gc_percent)}
                />
                <KV
                  label="DNA length"
                  value={`${job.validation_result.dna_length?.toLocaleString() ?? "—"} bp`}
                />
                <KV
                  label="ORFs detected"
                  value={job.validation_result.orfs?.length ?? 0}
                />
              </dl>
              {!!job.validation_result.warnings?.length && (
                <div className="mt-3 space-y-1">
                  {job.validation_result.warnings.map((w, i) => (
                    <p key={i} className="text-xs text-amber-700">
                      ⚠ {w}
                    </p>
                  ))}
                </div>
              )}
            </CardBody>
          </Card>
        )}

        {/* Protein & structure */}
        {job.analysis_result && (
          <Card>
            <CardHeader
              title="Protein & structure"
              icon={<Boxes className="h-4 w-4" />}
            />
            <CardBody>
              <dl className="divide-y divide-slate-100">
                <KV
                  label="BLAST status"
                  value={job.analysis_result.blast?.protein_status ?? "—"}
                />
                <KV
                  label="Top hit"
                  value={job.analysis_result.blast?.top_hit_name ?? "—"}
                />
                <KV
                  label="Identity"
                  value={formatPercent(
                    job.analysis_result.blast?.identity_percent
                  )}
                />
                <KV
                  label="UniProt"
                  value={
                    <span className="font-mono text-xs">
                      {job.analysis_result.annotation?.uniprot_id ?? "—"}
                    </span>
                  }
                />
                <KV
                  label="Foldable"
                  value={
                    <Badge
                      tone={
                        job.analysis_result.fold_check?.foldable
                          ? "green"
                          : "amber"
                      }
                    >
                      {job.analysis_result.fold_check?.foldable ? "Yes" : "No"}
                    </Badge>
                  }
                />
                <KV
                  label="Disordered"
                  value={formatPercent(
                    job.analysis_result.fold_check?.disordered_percent
                  )}
                />
                <KV
                  label="Structure source"
                  value={job.analysis_result.structure_3d?.source ?? "—"}
                />
                <KV
                  label="Mean pLDDT"
                  value={
                    job.analysis_result.structure_3d?.plddt_mean != null
                      ? job.analysis_result.structure_3d.plddt_mean.toFixed(1)
                      : "—"
                  }
                />
                <KV
                  label="Confidence"
                  value={
                    <Badge tone="blue">
                      {job.analysis_result.structure_3d?.confidence_level ?? "—"}
                    </Badge>
                  }
                />
              </dl>

              {!!job.analysis_result.annotation?.diseases?.length && (
                <div className="mt-3">
                  <p className="mb-1 text-xs font-medium text-slate-500">
                    Associated diseases
                  </p>
                  <div className="flex flex-wrap gap-1.5">
                    {job.analysis_result.annotation.diseases
                      .slice(0, 8)
                      .map((d, i) => (
                        <Badge key={i} tone="slate">
                          {d}
                        </Badge>
                      ))}
                  </div>
                </div>
              )}

              {/* 3D viewer placeholder */}
              <div className="mt-4 rounded-lg border border-dashed border-slate-300 bg-slate-50 p-4 text-center">
                <Layers className="mx-auto h-6 w-6 text-slate-400" />
                <p className="mt-2 text-sm font-medium text-slate-600">
                  3D structure viewer
                </p>
                <p className="mt-0.5 text-xs text-slate-400">
                  py3Dmol output served from{" "}
                  <span className="font-mono">
                    {job.analysis_result.structure_3d?.pdb_file_path || "—"}
                  </span>
                </p>
              </div>
            </CardBody>
          </Card>
        )}

        {/* Drug candidates */}
        {job.drug_result && (
          <Card className="lg:col-span-2">
            <CardHeader
              title="Drug candidates"
              description={
                job.drug_result.pocket
                  ? `Binding pocket score ${
                      job.drug_result.pocket.score ?? "—"
                    } · volume ${job.drug_result.pocket.volume ?? "—"} Å³`
                  : undefined
              }
              icon={<Pill className="h-4 w-4" />}
            />
            <Table>
              <THead>
                <TR>
                  <TH>ChEMBL ID</TH>
                  <TH>SMILES</TH>
                  <TH>Binding</TH>
                  <TH>Toxicity</TH>
                  <TH>Lipinski</TH>
                </TR>
              </THead>
              <TBody>
                {!job.drug_result.top_candidates?.length ? (
                  <EmptyRow colSpan={5} label="No candidates returned." />
                ) : (
                  job.drug_result.top_candidates.map((c, i) => (
                    <TR key={c.chembl_id || i}>
                      <TD className="font-mono text-xs font-medium text-slate-900">
                        {c.chembl_id || "—"}
                      </TD>
                      <TD>
                        <code
                          title={c.smiles}
                          className="font-mono text-xs text-slate-600"
                        >
                          {truncateMiddle(c.smiles || "", 36)}
                        </code>
                      </TD>
                      <TD className="tabular-nums">
                        {c.binding_score != null
                          ? c.binding_score.toFixed(2)
                          : "—"}
                      </TD>
                      <TD>
                        <Badge tone={toxicityTone(c.toxicity_risk)}>
                          {c.toxicity_risk || "—"}
                        </Badge>
                      </TD>
                      <TD>
                        <Badge tone={c.lipinski_pass ? "green" : "amber"}>
                          {c.lipinski_pass ? "Pass" : "Fail"}
                        </Badge>
                      </TD>
                    </TR>
                  ))
                )}
              </TBody>
            </Table>
            {job.drug_result.visualization_html_path && (
              <CardBody className="border-t border-slate-100">
                <p className="text-xs text-slate-500">
                  Interactive docking visualization served from{" "}
                  <span className="font-mono">
                    {job.drug_result.visualization_html_path}
                  </span>
                </p>
              </CardBody>
            )}
          </Card>
        )}

        {/* Histopathology */}
        {job.histopathology_result && (
          <Card className="lg:col-span-2">
            <CardHeader
              title="Histopathology"
              description={`Model: ${job.histopathology_result.model ?? "—"}`}
              icon={<Microscope className="h-4 w-4" />}
            />
            <CardBody className="space-y-4">
              {/* Track B is integrated end-to-end but its model is not trained yet. */}
              <div className="rounded-lg border border-amber-300 bg-amber-50 px-4 py-3 text-sm text-amber-800">
                <span className="font-semibold">Demo — model not trained yet.</span>{" "}
                The Track&nbsp;B histopathology model runs on untrained weights. These
                mutation predictions illustrate the pipeline only and are{" "}
                <span className="font-semibold">not scientifically valid</span>.
              </div>
              {!!job.histopathology_result.warnings?.length && (
                <InfoBanner>
                  {job.histopathology_result.warnings.map((w, i) => (
                    <p key={i} className="font-medium">
                      {w}
                    </p>
                  ))}
                </InfoBanner>
              )}

              <div className="grid grid-cols-1 gap-6 md:grid-cols-2">
                <div>
                  <dl className="divide-y divide-slate-100">
                    <KV
                      label="Patches kept"
                      value={`${
                        job.histopathology_result.patches_kept ?? "—"
                      } / ${job.histopathology_result.num_patches ?? "—"}`}
                    />
                    <KV
                      label="Mutations evaluated"
                      value={job.histopathology_result.mutations?.length ?? 0}
                    />
                  </dl>

                  <div className="mt-4 rounded-lg border border-dashed border-slate-300 bg-slate-50 p-4 text-center">
                    <Microscope className="mx-auto h-6 w-6 text-slate-400" />
                    <p className="mt-2 text-sm font-medium text-slate-600">
                      Grad-CAM overlay
                    </p>
                    <p className="mt-0.5 text-xs text-slate-400">
                      Attention heatmap served from{" "}
                      <span className="font-mono">
                        {basename(
                          job.histopathology_result.gradcam_overlay_path
                        ) || "—"}
                      </span>
                    </p>
                  </div>
                </div>

                <div>
                  <p className="mb-2 text-xs font-medium uppercase tracking-wide text-slate-500">
                    Per-gene mutation predictions
                  </p>
                  <div className="space-y-3">
                    {!job.histopathology_result.mutations?.length ? (
                      <p className="text-sm text-slate-500">
                        No predictions available.
                      </p>
                    ) : (
                      job.histopathology_result.mutations.map((m) => (
                        <div key={m.gene}>
                          <div className="mb-1 flex items-center justify-between">
                            <span className="flex items-center gap-2 text-sm font-medium text-slate-800">
                              <CircleDot
                                className={cn(
                                  "h-3.5 w-3.5",
                                  m.mutated
                                    ? "text-red-500"
                                    : "text-slate-300"
                                )}
                              />
                              {m.gene}
                            </span>
                            <Badge tone={m.mutated ? "red" : "slate"}>
                              {m.mutated ? "Mutated" : "Wild-type"}
                            </Badge>
                          </div>
                          <ConfidenceBar
                            value={m.confidence ?? m.probability}
                            label={`p = ${formatPercent(m.probability)}`}
                          />
                        </div>
                      ))
                    )}
                  </div>
                </div>
              </div>
            </CardBody>
          </Card>
        )}

        {/* Cross-modal fusion */}
        {job.fusion_result && (
          <Card className="lg:col-span-2">
            <CardHeader
              title="Cross-modal fusion"
              description="Concordance between genomic and visual signals."
              icon={<GitCompareArrows className="h-4 w-4" />}
              action={
                <div className="flex items-center gap-2">
                  <Badge tone="brand">
                    Agreement {formatPercent(job.fusion_result.overall_agreement)}
                  </Badge>
                  <Badge tone="slate">
                    κ ={" "}
                    {job.fusion_result.cohen_kappa != null
                      ? job.fusion_result.cohen_kappa.toFixed(2)
                      : "—"}
                  </Badge>
                </div>
              }
            />
            <Table>
              <THead>
                <TR>
                  <TH>Gene</TH>
                  <TH>Genomic</TH>
                  <TH>Visual</TH>
                  <TH>Agreement</TH>
                  <TH className="w-48">Combined confidence</TH>
                  <TH>Review</TH>
                </TR>
              </THead>
              <TBody>
                {!job.fusion_result.genes?.length ? (
                  <EmptyRow colSpan={6} label="No fusion results." />
                ) : (
                  job.fusion_result.genes.map((g) => (
                    <TR key={g.gene}>
                      <TD className="font-medium text-slate-900">{g.gene}</TD>
                      <TD className="text-slate-600">{g.genomic_signal}</TD>
                      <TD className="text-slate-600">
                        {g.visual_signal}
                        {g.visual_probability != null && (
                          <span className="ml-1 text-xs text-slate-400">
                            ({formatPercent(g.visual_probability)})
                          </span>
                        )}
                      </TD>
                      <TD>
                        <Badge tone={agreementTone(g.agreement)}>
                          {g.agreement}
                        </Badge>
                      </TD>
                      <TD>
                        <ConfidenceBar value={g.combined_confidence} />
                      </TD>
                      <TD>
                        {g.flag_for_review ? (
                          <Badge tone="red">Flagged</Badge>
                        ) : (
                          <Badge tone="green">OK</Badge>
                        )}
                      </TD>
                    </TR>
                  ))
                )}
              </TBody>
            </Table>
            {!!job.fusion_result.flagged_genes?.length && (
              <CardBody className="border-t border-slate-100">
                <p className="text-xs text-slate-500">
                  Flagged for review:{" "}
                  <span className="font-medium text-slate-700">
                    {job.fusion_result.flagged_genes.join(", ")}
                  </span>
                </p>
              </CardBody>
            )}
          </Card>
        )}
      </div>

      {/* When completed but no stage results returned */}
      {isCompleted &&
        !job.validation_result &&
        !job.analysis_result &&
        !job.drug_result &&
        !job.histopathology_result &&
        !job.fusion_result && (
          <Card className="mt-6">
            <CardBody className="py-10 text-center text-sm text-slate-500">
              This analysis completed but returned no stage results.
            </CardBody>
          </Card>
        )}

      {isCompleted && (
        <div className="mt-6 flex justify-end">
          <Button onClick={onGenerateReport} loading={generating}>
            <FileText className="h-4 w-4" /> Generate report
            <ExternalLink className="h-3.5 w-3.5 opacity-70" />
          </Button>
        </div>
      )}

      <Modal
        open={confirmDelete}
        onClose={() => setConfirmDelete(false)}
        title="Delete analysis"
        description="This action cannot be undone."
        footer={
          <>
            <Button variant="outline" onClick={() => setConfirmDelete(false)}>
              Cancel
            </Button>
            <Button variant="danger" onClick={onDelete} loading={deleting}>
              Delete analysis
            </Button>
          </>
        }
      >
        <p className="text-sm text-slate-600">
          Are you sure you want to delete <strong>{job.name}</strong> and all of
          its results?
        </p>
      </Modal>
    </div>
  );
}
