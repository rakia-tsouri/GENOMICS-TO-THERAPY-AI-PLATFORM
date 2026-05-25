"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import {
  FlaskConical,
  CheckCircle2,
  FolderKanban,
  FileText,
  Plus,
  ArrowRight,
} from "lucide-react";
import { get, ApiError } from "@/lib/api";
import type { JobSummary, Project, Report } from "@/lib/types";
import { useAuth } from "@/lib/auth";
import PageHeader from "@/components/PageHeader";
import StatChip from "@/components/StatChip";
import StatusBadge from "@/components/StatusBadge";
import Button from "@/components/Button";
import Card, { CardHeader, CardBody } from "@/components/Card";
import { Table, THead, TBody, TR, TH, TD, EmptyRow } from "@/components/Table";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import { formatDate } from "@/lib/utils";

export default function DashboardPage() {
  const { user } = useAuth();
  const [jobs, setJobs] = useState<JobSummary[]>([]);
  const [projects, setProjects] = useState<Project[]>([]);
  const [reports, setReports] = useState<Report[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    async function load() {
      setLoading(true);
      setError(null);
      try {
        const [j, p, r] = await Promise.all([
          get<JobSummary[]>("/jobs"),
          get<Project[]>("/projects"),
          get<Report[]>("/reports"),
        ]);
        if (cancelled) return;
        setJobs(j);
        setProjects(p);
        setReports(r);
      } catch (err) {
        if (cancelled) return;
        setError(
          err instanceof ApiError ? err.message : "Failed to load dashboard."
        );
      } finally {
        if (!cancelled) setLoading(false);
      }
    }
    void load();
    return () => {
      cancelled = true;
    };
  }, []);

  const completed = jobs.filter((j) => j.status === "completed").length;
  const recent = [...jobs]
    .sort(
      (a, b) =>
        new Date(b.created_at).getTime() - new Date(a.created_at).getTime()
    )
    .slice(0, 6);

  if (loading) return <PageSpinner label="Loading dashboard…" />;

  return (
    <div>
      <PageHeader
        title={`Welcome${user?.full_name ? `, ${user.full_name.split(" ")[0]}` : ""}`}
        description="Overview of your analyses, projects and reports."
        action={
          <Link href="/jobs/new">
            <Button>
              <Plus className="h-4 w-4" />
              New analysis
            </Button>
          </Link>
        }
      />

      <ErrorBanner message={error} className="mb-6" />

      <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">
        <StatChip
          label="Total analyses"
          value={jobs.length}
          icon={<FlaskConical className="h-4 w-4" />}
        />
        <StatChip
          label="Completed"
          value={completed}
          icon={<CheckCircle2 className="h-4 w-4" />}
          hint={
            jobs.length > 0
              ? `${Math.round((completed / jobs.length) * 100)}% success rate`
              : undefined
          }
        />
        <StatChip
          label="Projects"
          value={projects.length}
          icon={<FolderKanban className="h-4 w-4" />}
        />
        <StatChip
          label="Reports"
          value={reports.length}
          icon={<FileText className="h-4 w-4" />}
        />
      </div>

      <div className="mt-6">
        <Card>
          <CardHeader
            title="Recent analyses"
            description="Your most recently submitted pipeline runs."
            icon={<FlaskConical className="h-4 w-4" />}
            action={
              <Link
                href="/jobs"
                className="inline-flex items-center gap-1 text-sm font-medium text-brand-600 hover:text-brand-700"
              >
                View all <ArrowRight className="h-4 w-4" />
              </Link>
            }
          />
          <Table>
            <THead>
              <TR>
                <TH>Name</TH>
                <TH>Status</TH>
                <TH>Step</TH>
                <TH>Gene</TH>
                <TH>Created</TH>
                <TH className="text-right">Actions</TH>
              </TR>
            </THead>
            <TBody>
              {recent.length === 0 ? (
                <EmptyRow colSpan={6} label="No analyses yet. Start your first run." />
              ) : (
                recent.map((job) => (
                  <TR key={job.id}>
                    <TD className="font-medium text-slate-900">
                      <Link href={`/jobs/${job.id}`} className="hover:text-brand-600">
                        {job.name}
                      </Link>
                    </TD>
                    <TD>
                      <StatusBadge status={job.status} />
                    </TD>
                    <TD className="text-slate-500">
                      {job.current_step || "—"}
                    </TD>
                    <TD className="font-mono text-xs">{job.gene_id || "—"}</TD>
                    <TD className="text-slate-500">{formatDate(job.created_at)}</TD>
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
        </Card>
      </div>
    </div>
  );
}
