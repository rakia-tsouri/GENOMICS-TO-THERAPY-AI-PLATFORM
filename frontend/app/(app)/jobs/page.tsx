"use client";

import { useEffect, useMemo, useState } from "react";
import Link from "next/link";
import { FlaskConical, Plus } from "lucide-react";
import { get, ApiError } from "@/lib/api";
import type { JobStatus, JobSummary, Project } from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Card, { CardBody } from "@/components/Card";
import StatusBadge from "@/components/StatusBadge";
import { Table, THead, TBody, TR, TH, TD, EmptyRow } from "@/components/Table";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import { cn, formatDate } from "@/lib/utils";

const FILTERS: { key: "all" | JobStatus; label: string }[] = [
  { key: "all", label: "All" },
  { key: "pending", label: "Pending" },
  { key: "running", label: "Running" },
  { key: "completed", label: "Completed" },
  { key: "failed", label: "Failed" },
];

export default function JobsPage() {
  const [jobs, setJobs] = useState<JobSummary[]>([]);
  const [projects, setProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [filter, setFilter] = useState<"all" | JobStatus>("all");

  useEffect(() => {
    let cancelled = false;
    async function load() {
      setLoading(true);
      setError(null);
      try {
        const [j, p] = await Promise.all([
          get<JobSummary[]>("/jobs"),
          get<Project[]>("/projects"),
        ]);
        if (cancelled) return;
        setJobs(j);
        setProjects(p);
      } catch (err) {
        if (cancelled) return;
        setError(err instanceof ApiError ? err.message : "Failed to load analyses.");
      } finally {
        if (!cancelled) setLoading(false);
      }
    }
    void load();
    return () => {
      cancelled = true;
    };
  }, []);

  const projectName = useMemo(() => {
    const map = new Map<number, string>();
    projects.forEach((p) => map.set(p.id, p.name));
    return map;
  }, [projects]);

  const filtered = useMemo(() => {
    const list = filter === "all" ? jobs : jobs.filter((j) => j.status === filter);
    return [...list].sort(
      (a, b) =>
        new Date(b.created_at).getTime() - new Date(a.created_at).getTime()
    );
  }, [jobs, filter]);

  return (
    <div>
      <PageHeader
        title="Analyses"
        description="All pipeline runs across your projects."
        action={
          <Link href="/jobs/new">
            <Button>
              <Plus className="h-4 w-4" /> New analysis
            </Button>
          </Link>
        }
      />

      <ErrorBanner message={error} className="mb-6" />

      <div className="mb-4 flex flex-wrap gap-2">
        {FILTERS.map((f) => {
          const count =
            f.key === "all"
              ? jobs.length
              : jobs.filter((j) => j.status === f.key).length;
          return (
            <button
              key={f.key}
              onClick={() => setFilter(f.key)}
              className={cn(
                "rounded-full border px-3 py-1.5 text-sm font-medium transition-colors",
                filter === f.key
                  ? "border-brand-300 bg-brand-50 text-brand-700"
                  : "border-slate-200 bg-white text-slate-600 hover:bg-slate-50"
              )}
            >
              {f.label}
              <span className="ml-1.5 text-xs text-slate-400">{count}</span>
            </button>
          );
        })}
      </div>

      {loading ? (
        <PageSpinner label="Loading analyses…" />
      ) : (
        <Card>
          {filtered.length === 0 && filter === "all" && jobs.length === 0 ? (
            <CardBody className="flex flex-col items-center gap-3 py-14 text-center">
              <span className="flex h-12 w-12 items-center justify-center rounded-full bg-brand-50 text-brand-600">
                <FlaskConical className="h-6 w-6" />
              </span>
              <div>
                <p className="font-medium text-slate-900">No analyses yet</p>
                <p className="text-sm text-slate-500">
                  Launch your first end-to-end pipeline run.
                </p>
              </div>
              <Link href="/jobs/new">
                <Button>
                  <Plus className="h-4 w-4" /> New analysis
                </Button>
              </Link>
            </CardBody>
          ) : (
            <Table>
              <THead>
                <TR>
                  <TH>Name</TH>
                  <TH>Status</TH>
                  <TH>Step</TH>
                  <TH>Project</TH>
                  <TH>Gene</TH>
                  <TH>Created</TH>
                  <TH className="text-right">Actions</TH>
                </TR>
              </THead>
              <TBody>
                {filtered.length === 0 ? (
                  <EmptyRow colSpan={7} label="No analyses match this filter." />
                ) : (
                  filtered.map((job) => (
                    <TR key={job.id}>
                      <TD className="font-medium text-slate-900">
                        <Link
                          href={`/jobs/${job.id}`}
                          className="hover:text-brand-600"
                        >
                          {job.name}
                        </Link>
                      </TD>
                      <TD>
                        <StatusBadge status={job.status} />
                      </TD>
                      <TD className="text-slate-500">{job.current_step || "—"}</TD>
                      <TD className="text-slate-500">
                        {job.project_id
                          ? projectName.get(job.project_id) ??
                            `#${job.project_id}`
                          : "—"}
                      </TD>
                      <TD className="font-mono text-xs">{job.gene_id || "—"}</TD>
                      <TD className="text-slate-500">
                        {formatDate(job.created_at)}
                      </TD>
                      <TD className="text-right">
                        <Link
                          href={`/jobs/${job.id}`}
                          className="text-sm font-medium text-brand-600 hover:text-brand-700"
                        >
                          Open
                        </Link>
                      </TD>
                    </TR>
                  ))
                )}
              </TBody>
            </Table>
          )}
        </Card>
      )}
    </div>
  );
}
