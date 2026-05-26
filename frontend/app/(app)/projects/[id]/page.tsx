"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { useParams, useRouter } from "next/navigation";
import { ArrowLeft, Plus, Trash2, Pencil, FlaskConical } from "lucide-react";
import { get, patch, del, ApiError } from "@/lib/api";
import type { JobSummary, Project } from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Card, { CardHeader } from "@/components/Card";
import Badge from "@/components/Badge";
import Modal from "@/components/Modal";
import { Field, Input, Textarea } from "@/components/Input";
import StatusBadge from "@/components/StatusBadge";
import { Table, THead, TBody, TR, TH, TD, EmptyRow } from "@/components/Table";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import { formatDate } from "@/lib/utils";

export default function ProjectDetailPage() {
  const params = useParams<{ id: string }>();
  const router = useRouter();
  const projectId = Number(params.id);

  const [project, setProject] = useState<Project | null>(null);
  const [jobs, setJobs] = useState<JobSummary[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const [editOpen, setEditOpen] = useState(false);
  const [editForm, setEditForm] = useState({
    name: "",
    description: "",
    cancer_type: "",
  });
  const [saving, setSaving] = useState(false);
  const [deleting, setDeleting] = useState(false);
  const [confirmDelete, setConfirmDelete] = useState(false);

  useEffect(() => {
    let cancelled = false;
    async function load() {
      setLoading(true);
      setError(null);
      try {
        const [p, j] = await Promise.all([
          get<Project>(`/projects/${projectId}`),
          get<JobSummary[]>(`/jobs?project_id=${projectId}`),
        ]);
        if (cancelled) return;
        setProject(p);
        setJobs(j);
        setEditForm({
          name: p.name,
          description: p.description,
          cancer_type: p.cancer_type,
        });
      } catch (err) {
        if (cancelled) return;
        setError(
          err instanceof ApiError ? err.message : "Failed to load project."
        );
      } finally {
        if (!cancelled) setLoading(false);
      }
    }
    if (!Number.isNaN(projectId)) void load();
    return () => {
      cancelled = true;
    };
  }, [projectId]);

  async function onSave(e: React.FormEvent) {
    e.preventDefault();
    setSaving(true);
    setError(null);
    try {
      const updated = await patch<Project>(`/projects/${projectId}`, editForm);
      setProject(updated);
      setEditOpen(false);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to save changes.");
    } finally {
      setSaving(false);
    }
  }

  async function onDelete() {
    setDeleting(true);
    setError(null);
    try {
      await del(`/projects/${projectId}`);
      router.replace("/projects");
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to delete project.");
      setDeleting(false);
    }
  }

  if (loading) return <PageSpinner label="Loading project…" />;

  if (!project) {
    return (
      <div>
        <ErrorBanner message={error || "Project not found."} />
        <Link
          href="/projects"
          className="mt-4 inline-flex items-center gap-1 text-sm text-brand-600"
        >
          <ArrowLeft className="h-4 w-4" /> Back to projects
        </Link>
      </div>
    );
  }

  return (
    <div>
      <Link
        href="/projects"
        className="mb-4 inline-flex items-center gap-1 text-sm text-slate-500 hover:text-slate-700"
      >
        <ArrowLeft className="h-4 w-4" /> Projects
      </Link>

      <PageHeader
        title={
          <span className="flex items-center gap-3">
            {project.name}
            <Badge tone="brand">{project.cancer_type}</Badge>
          </span>
        }
        description={project.description || "No description provided."}
        action={
          <div className="flex gap-2">
            <Button variant="outline" onClick={() => setEditOpen(true)}>
              <Pencil className="h-4 w-4" /> Edit
            </Button>
            <Button variant="danger" onClick={() => setConfirmDelete(true)}>
              <Trash2 className="h-4 w-4" /> Delete
            </Button>
          </div>
        }
      />

      <ErrorBanner message={error} className="mb-6" />

      <Card>
        <CardHeader
          title="Analyses in this project"
          icon={<FlaskConical className="h-4 w-4" />}
          action={
            <Link href={`/jobs/new?project_id=${projectId}`}>
              <Button size="sm">
                <Plus className="h-4 w-4" /> New analysis
              </Button>
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
            {jobs.length === 0 ? (
              <EmptyRow colSpan={6} label="No analyses in this project yet." />
            ) : (
              jobs.map((job) => (
                <TR key={job.id}>
                  <TD className="font-medium text-slate-900">
                    <Link href={`/jobs/${job.id}`} className="hover:text-brand-600">
                      {job.name}
                    </Link>
                  </TD>
                  <TD>
                    <StatusBadge status={job.status} />
                  </TD>
                  <TD className="text-slate-500">{job.current_step || "—"}</TD>
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

      <Modal
        open={editOpen}
        onClose={() => setEditOpen(false)}
        title="Edit project"
        footer={
          <>
            <Button variant="outline" onClick={() => setEditOpen(false)}>
              Cancel
            </Button>
            <Button form="edit-project-form" type="submit" loading={saving}>
              Save changes
            </Button>
          </>
        }
      >
        <form id="edit-project-form" onSubmit={onSave} className="space-y-4">
          <Field label="Name" htmlFor="e-name" required>
            <Input
              id="e-name"
              required
              value={editForm.name}
              onChange={(e) => setEditForm({ ...editForm, name: e.target.value })}
            />
          </Field>
          <Field label="Cancer type" htmlFor="e-type">
            <Input
              id="e-type"
              value={editForm.cancer_type}
              onChange={(e) =>
                setEditForm({ ...editForm, cancer_type: e.target.value })
              }
            />
          </Field>
          <Field label="Description" htmlFor="e-desc">
            <Textarea
              id="e-desc"
              className="min-h-[90px] font-sans text-sm"
              value={editForm.description}
              onChange={(e) =>
                setEditForm({ ...editForm, description: e.target.value })
              }
            />
          </Field>
        </form>
      </Modal>

      <Modal
        open={confirmDelete}
        onClose={() => setConfirmDelete(false)}
        title="Delete project"
        description="This action cannot be undone."
        footer={
          <>
            <Button variant="outline" onClick={() => setConfirmDelete(false)}>
              Cancel
            </Button>
            <Button variant="danger" onClick={onDelete} loading={deleting}>
              Delete project
            </Button>
          </>
        }
      >
        <p className="text-sm text-slate-600">
          Are you sure you want to delete <strong>{project.name}</strong>? Its
          analyses will not be re-assigned.
        </p>
      </Modal>
    </div>
  );
}
