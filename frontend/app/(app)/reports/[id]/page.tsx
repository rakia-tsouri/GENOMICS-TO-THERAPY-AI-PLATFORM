"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { useParams } from "next/navigation";
import { ArrowLeft, Download, FlaskConical } from "lucide-react";
import { get, fetchBlob, ApiError } from "@/lib/api";
import type { Report } from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Card, { CardHeader, CardBody } from "@/components/Card";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import { formatDate } from "@/lib/utils";

/** Recursively render an arbitrary JSON summary in a readable way. */
function SummaryValue({ value }: { value: unknown }) {
  if (value === null || value === undefined) {
    return <span className="text-slate-400">—</span>;
  }
  if (typeof value === "boolean") {
    return <span className="font-medium">{value ? "Yes" : "No"}</span>;
  }
  if (typeof value === "number" || typeof value === "string") {
    return <span className="break-words">{String(value)}</span>;
  }
  if (Array.isArray(value)) {
    if (value.length === 0) return <span className="text-slate-400">[]</span>;
    const primitives = value.every(
      (v) => typeof v === "string" || typeof v === "number" || typeof v === "boolean"
    );
    if (primitives) {
      return (
        <div className="flex flex-wrap gap-1.5">
          {value.map((v, i) => (
            <span
              key={i}
              className="rounded-md bg-slate-100 px-2 py-0.5 text-xs text-slate-700"
            >
              {String(v)}
            </span>
          ))}
        </div>
      );
    }
    return (
      <div className="space-y-2">
        {value.map((v, i) => (
          <div
            key={i}
            className="rounded-lg border border-slate-100 bg-slate-50/60 p-3"
          >
            <SummaryValue value={v} />
          </div>
        ))}
      </div>
    );
  }
  if (typeof value === "object") {
    return <SummaryObject obj={value as Record<string, unknown>} />;
  }
  return <span>{String(value)}</span>;
}

function SummaryObject({ obj }: { obj: Record<string, unknown> }) {
  const entries = Object.entries(obj);
  if (entries.length === 0) return <span className="text-slate-400">{"{}"}</span>;
  return (
    <dl className="space-y-2">
      {entries.map(([key, val]) => (
        <div
          key={key}
          className="grid grid-cols-1 gap-1 sm:grid-cols-3 sm:gap-3"
        >
          <dt className="text-xs font-medium uppercase tracking-wide text-slate-500 sm:pt-0.5">
            {key.replace(/_/g, " ")}
          </dt>
          <dd className="text-sm text-slate-800 sm:col-span-2">
            <SummaryValue value={val} />
          </dd>
        </div>
      ))}
    </dl>
  );
}

export default function ReportDetailPage() {
  const params = useParams<{ id: string }>();
  const reportId = Number(params.id);

  const [report, setReport] = useState<Report | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [exporting, setExporting] = useState(false);

  useEffect(() => {
    let cancelled = false;
    async function load() {
      setLoading(true);
      setError(null);
      try {
        const data = await get<Report>(`/reports/${reportId}`);
        if (!cancelled) setReport(data);
      } catch (err) {
        if (!cancelled)
          setError(
            err instanceof ApiError ? err.message : "Failed to load report."
          );
      } finally {
        if (!cancelled) setLoading(false);
      }
    }
    if (!Number.isNaN(reportId)) void load();
    return () => {
      cancelled = true;
    };
  }, [reportId]);

  async function onExport() {
    if (!report) return;
    setExporting(true);
    setError(null);
    try {
      const { blob } = await fetchBlob(`/reports/${report.id}/export`);
      const url = URL.createObjectURL(blob);
      window.open(url, "_blank", "noopener,noreferrer");
      setTimeout(() => URL.revokeObjectURL(url), 60_000);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to export report.");
    } finally {
      setExporting(false);
    }
  }

  if (loading) return <PageSpinner label="Loading report…" />;

  if (!report) {
    return (
      <div>
        <ErrorBanner message={error || "Report not found."} />
        <Link
          href="/reports"
          className="mt-4 inline-flex items-center gap-1 text-sm text-brand-600"
        >
          <ArrowLeft className="h-4 w-4" /> Back to reports
        </Link>
      </div>
    );
  }

  return (
    <div>
      <Link
        href="/reports"
        className="mb-4 inline-flex items-center gap-1 text-sm text-slate-500 hover:text-slate-700"
      >
        <ArrowLeft className="h-4 w-4" /> Reports
      </Link>

      <PageHeader
        title={report.title || `Report #${report.id}`}
        description={`Generated ${formatDate(report.created_at)}`}
        action={
          <Button onClick={onExport} loading={exporting}>
            <Download className="h-4 w-4" /> Export HTML
          </Button>
        }
      />

      <ErrorBanner message={error} className="mb-6" />

      <div className="mb-6">
        <Link
          href={`/jobs/${report.job_id}`}
          className="inline-flex items-center gap-1.5 rounded-lg border border-slate-200 bg-white px-3 py-2 text-sm font-medium text-brand-600 hover:bg-slate-50"
        >
          <FlaskConical className="h-4 w-4" /> View source analysis (Job #
          {report.job_id})
        </Link>
      </div>

      <Card>
        <CardHeader
          title="Report summary"
          description="Structured summary captured at generation time."
        />
        <CardBody>
          {report.summary && Object.keys(report.summary).length > 0 ? (
            <SummaryObject obj={report.summary} />
          ) : (
            <p className="text-sm text-slate-500">
              No structured summary available. Use{" "}
              <strong>Export HTML</strong> for the full report document.
            </p>
          )}
        </CardBody>
      </Card>
    </div>
  );
}
