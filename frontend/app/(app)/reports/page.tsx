"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { FileText, Download, Trash2 } from "lucide-react";
import { get, del, fetchBlob, ApiError, API_BASE, API_PREFIX } from "@/lib/api";
import type { Report } from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Card, { CardBody } from "@/components/Card";
import { Table, THead, TBody, TR, TH, TD, EmptyRow } from "@/components/Table";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import { formatDate } from "@/lib/utils";

export default function ReportsPage() {
  const [reports, setReports] = useState<Report[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [exportingId, setExportingId] = useState<number | null>(null);
  const [deletingId, setDeletingId] = useState<number | null>(null);

  async function load() {
    setLoading(true);
    setError(null);
    try {
      const data = await get<Report[]>("/reports");
      setReports(data);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to load reports.");
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    void load();
  }, []);

  async function onExport(report: Report) {
    setExportingId(report.id);
    setError(null);
    try {
      // Export needs the bearer token, so we fetch as a blob and open it.
      const { blob } = await fetchBlob(`/reports/${report.id}/export`);
      const url = URL.createObjectURL(blob);
      window.open(url, "_blank", "noopener,noreferrer");
      // Revoke after a delay so the new tab can load it.
      setTimeout(() => URL.revokeObjectURL(url), 60_000);
    } catch (err) {
      setError(
        err instanceof ApiError ? err.message : "Failed to export report."
      );
    } finally {
      setExportingId(null);
    }
  }

  async function onDelete(report: Report) {
    setDeletingId(report.id);
    setError(null);
    try {
      await del(`/reports/${report.id}`);
      setReports((prev) => prev.filter((r) => r.id !== report.id));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to delete report.");
    } finally {
      setDeletingId(null);
    }
  }

  return (
    <div>
      <PageHeader
        title="Reports"
        description="Shareable summaries generated from completed analyses."
      />

      <ErrorBanner message={error} className="mb-6" />

      {loading ? (
        <PageSpinner label="Loading reports…" />
      ) : (
        <Card>
          {reports.length === 0 ? (
            <CardBody className="flex flex-col items-center gap-3 py-14 text-center">
              <span className="flex h-12 w-12 items-center justify-center rounded-full bg-brand-50 text-brand-600">
                <FileText className="h-6 w-6" />
              </span>
              <div>
                <p className="font-medium text-slate-900">No reports yet</p>
                <p className="text-sm text-slate-500">
                  Generate a report from any completed analysis.
                </p>
              </div>
              <Link href="/jobs">
                <Button variant="outline">Go to analyses</Button>
              </Link>
            </CardBody>
          ) : (
            <Table>
              <THead>
                <TR>
                  <TH>Title</TH>
                  <TH>Analysis</TH>
                  <TH>Created</TH>
                  <TH className="text-right">Actions</TH>
                </TR>
              </THead>
              <TBody>
                {reports.map((r) => (
                  <TR key={r.id}>
                    <TD className="font-medium text-slate-900">
                      <Link
                        href={`/reports/${r.id}`}
                        className="hover:text-brand-600"
                      >
                        {r.title || `Report #${r.id}`}
                      </Link>
                    </TD>
                    <TD>
                      <Link
                        href={`/jobs/${r.job_id}`}
                        className="text-sm text-brand-600 hover:text-brand-700"
                      >
                        Job #{r.job_id}
                      </Link>
                    </TD>
                    <TD className="text-slate-500">{formatDate(r.created_at)}</TD>
                    <TD>
                      <div className="flex items-center justify-end gap-2">
                        <Button
                          size="sm"
                          variant="outline"
                          loading={exportingId === r.id}
                          onClick={() => onExport(r)}
                        >
                          <Download className="h-4 w-4" /> Export HTML
                        </Button>
                        <Button
                          size="sm"
                          variant="ghost"
                          loading={deletingId === r.id}
                          onClick={() => onDelete(r)}
                          aria-label="Delete report"
                        >
                          <Trash2 className="h-4 w-4 text-red-500" />
                        </Button>
                      </div>
                    </TD>
                  </TR>
                ))}
              </TBody>
            </Table>
          )}
        </Card>
      )}

      <p className="mt-4 text-xs text-slate-400">
        Export opens the authenticated HTML document (
        <span className="font-mono">
          {API_BASE}
          {API_PREFIX}/reports/&#123;id&#125;/export
        </span>
        ) in a new tab.
      </p>
    </div>
  );
}
