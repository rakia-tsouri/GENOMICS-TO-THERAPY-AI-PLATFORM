"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { FolderKanban, Plus } from "lucide-react";
import { get, post, ApiError } from "@/lib/api";
import type { CreateProjectInput, Project } from "@/lib/types";
import PageHeader from "@/components/PageHeader";
import Button from "@/components/Button";
import Modal from "@/components/Modal";
import { Field, Input, Textarea } from "@/components/Input";
import Card, { CardBody } from "@/components/Card";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import Badge from "@/components/Badge";
import { formatDate } from "@/lib/utils";

const CANCER_TYPES = [
  "Glioma",
  "Breast",
  "Lung (NSCLC)",
  "Colorectal",
  "Melanoma",
  "Pancreatic",
  "Leukemia (AML)",
  "Prostate",
  "Other",
];

export default function ProjectsPage() {
  const [projects, setProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const [open, setOpen] = useState(false);
  const [form, setForm] = useState<CreateProjectInput>({
    name: "",
    description: "",
    cancer_type: CANCER_TYPES[0],
  });
  const [submitting, setSubmitting] = useState(false);
  const [formError, setFormError] = useState<string | null>(null);

  async function load() {
    setLoading(true);
    setError(null);
    try {
      const data = await get<Project[]>("/projects");
      setProjects(data);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to load projects.");
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    void load();
  }, []);

  async function onCreate(e: React.FormEvent) {
    e.preventDefault();
    setFormError(null);
    setSubmitting(true);
    try {
      const created = await post<Project>("/projects", form);
      setProjects((prev) => [created, ...prev]);
      setOpen(false);
      setForm({ name: "", description: "", cancer_type: CANCER_TYPES[0] });
    } catch (err) {
      setFormError(
        err instanceof ApiError ? err.message : "Failed to create project."
      );
    } finally {
      setSubmitting(false);
    }
  }

  return (
    <div>
      <PageHeader
        title="Projects"
        description="Group related analyses under a project."
        action={
          <Button onClick={() => setOpen(true)}>
            <Plus className="h-4 w-4" />
            New project
          </Button>
        }
      />

      <ErrorBanner message={error} className="mb-6" />

      {loading ? (
        <PageSpinner label="Loading projects…" />
      ) : projects.length === 0 ? (
        <Card>
          <CardBody className="flex flex-col items-center gap-3 py-14 text-center">
            <span className="flex h-12 w-12 items-center justify-center rounded-full bg-brand-50 text-brand-600">
              <FolderKanban className="h-6 w-6" />
            </span>
            <div>
              <p className="font-medium text-slate-900">No projects yet</p>
              <p className="text-sm text-slate-500">
                Create a project to organize your analyses.
              </p>
            </div>
            <Button onClick={() => setOpen(true)}>
              <Plus className="h-4 w-4" />
              New project
            </Button>
          </CardBody>
        </Card>
      ) : (
        <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3">
          {projects.map((p) => (
            <Link key={p.id} href={`/projects/${p.id}`}>
              <Card className="h-full transition-shadow hover:shadow-md">
                <CardBody>
                  <div className="flex items-start justify-between gap-2">
                    <h3 className="font-semibold text-slate-900">{p.name}</h3>
                    <Badge tone="brand">{p.cancer_type}</Badge>
                  </div>
                  <p className="mt-2 line-clamp-3 text-sm text-slate-500">
                    {p.description || "No description provided."}
                  </p>
                  <p className="mt-4 text-xs text-slate-400">
                    Created {formatDate(p.created_at)}
                  </p>
                </CardBody>
              </Card>
            </Link>
          ))}
        </div>
      )}

      <Modal
        open={open}
        onClose={() => setOpen(false)}
        title="New project"
        description="Projects help you organize related analyses."
        footer={
          <>
            <Button variant="outline" onClick={() => setOpen(false)}>
              Cancel
            </Button>
            <Button form="new-project-form" type="submit" loading={submitting}>
              Create project
            </Button>
          </>
        }
      >
        <form id="new-project-form" onSubmit={onCreate} className="space-y-4">
          <ErrorBanner message={formError} />
          <Field label="Name" htmlFor="p-name" required>
            <Input
              id="p-name"
              required
              value={form.name}
              onChange={(e) => setForm({ ...form, name: e.target.value })}
              placeholder="e.g. TCGA-GBM cohort"
            />
          </Field>
          <Field label="Cancer type" htmlFor="p-type" required>
            <select
              id="p-type"
              className="w-full rounded-lg border border-slate-300 bg-white px-3 py-2 text-sm shadow-sm focus:border-brand-500 focus:outline-none focus:ring-2 focus:ring-brand-200"
              value={form.cancer_type}
              onChange={(e) => setForm({ ...form, cancer_type: e.target.value })}
            >
              {CANCER_TYPES.map((t) => (
                <option key={t} value={t}>
                  {t}
                </option>
              ))}
            </select>
          </Field>
          <Field label="Description" htmlFor="p-desc">
            <Textarea
              id="p-desc"
              className="min-h-[90px] font-sans text-sm"
              value={form.description}
              onChange={(e) => setForm({ ...form, description: e.target.value })}
              placeholder="What is this project about?"
            />
          </Field>
        </form>
      </Modal>
    </div>
  );
}
