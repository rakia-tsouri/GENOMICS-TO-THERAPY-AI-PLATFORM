"use client";

import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { ShieldAlert, Trash2 } from "lucide-react";
import { get, patch, del, ApiError } from "@/lib/api";
import type { Role, User } from "@/lib/types";
import { useAuth } from "@/lib/auth";
import PageHeader from "@/components/PageHeader";
import Card from "@/components/Card";
import Button from "@/components/Button";
import Badge from "@/components/Badge";
import Modal from "@/components/Modal";
import { Select } from "@/components/Input";
import { Table, THead, TBody, TR, TH, TD, EmptyRow } from "@/components/Table";
import { ErrorBanner, PageSpinner } from "@/components/Feedback";
import { formatDate } from "@/lib/utils";

export default function AdminUsersPage() {
  const { user: me, isAdmin, loading: authLoading } = useAuth();
  const router = useRouter();

  const [users, setUsers] = useState<User[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [busyId, setBusyId] = useState<number | null>(null);
  const [toDelete, setToDelete] = useState<User | null>(null);
  const [deleting, setDeleting] = useState(false);

  useEffect(() => {
    if (!authLoading && !isAdmin) {
      router.replace("/dashboard");
    }
  }, [authLoading, isAdmin, router]);

  async function load() {
    setLoading(true);
    setError(null);
    try {
      const data = await get<User[]>("/users");
      setUsers(data);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to load users.");
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    if (isAdmin) void load();
  }, [isAdmin]);

  async function updateUser(id: number, changes: Partial<User>) {
    setBusyId(id);
    setError(null);
    try {
      const updated = await patch<User>(`/users/${id}`, changes);
      setUsers((prev) => prev.map((u) => (u.id === id ? updated : u)));
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to update user.");
    } finally {
      setBusyId(null);
    }
  }

  async function onDelete() {
    if (!toDelete) return;
    setDeleting(true);
    setError(null);
    try {
      await del(`/users/${toDelete.id}`);
      setUsers((prev) => prev.filter((u) => u.id !== toDelete.id));
      setToDelete(null);
    } catch (err) {
      setError(err instanceof ApiError ? err.message : "Failed to delete user.");
    } finally {
      setDeleting(false);
    }
  }

  if (authLoading || !isAdmin) return <PageSpinner />;

  return (
    <div>
      <PageHeader
        title="User management"
        description="Manage roles and access for everyone on the platform."
      />

      <ErrorBanner message={error} className="mb-6" />

      {loading ? (
        <PageSpinner label="Loading users…" />
      ) : (
        <Card>
          <Table>
            <THead>
              <TR>
                <TH>User</TH>
                <TH>Role</TH>
                <TH>Status</TH>
                <TH>Joined</TH>
                <TH className="text-right">Actions</TH>
              </TR>
            </THead>
            <TBody>
              {users.length === 0 ? (
                <EmptyRow colSpan={5} label="No users found." />
              ) : (
                users.map((u) => {
                  const isSelf = me?.id === u.id;
                  const busy = busyId === u.id;
                  return (
                    <TR key={u.id}>
                      <TD>
                        <div>
                          <p className="font-medium text-slate-900">
                            {u.full_name || "—"}
                            {isSelf && (
                              <span className="ml-2 text-xs text-slate-400">
                                (you)
                              </span>
                            )}
                          </p>
                          <p className="text-xs text-slate-500">{u.email}</p>
                        </div>
                      </TD>
                      <TD>
                        <Select
                          value={u.role}
                          disabled={busy || isSelf}
                          className="h-8 w-36 py-1 text-xs"
                          onChange={(e) =>
                            updateUser(u.id, { role: e.target.value as Role })
                          }
                        >
                          <option value="researcher">researcher</option>
                          <option value="admin">admin</option>
                        </Select>
                      </TD>
                      <TD>
                        <button
                          type="button"
                          disabled={busy || isSelf}
                          onClick={() =>
                            updateUser(u.id, { is_active: !u.is_active })
                          }
                          className="disabled:cursor-not-allowed disabled:opacity-60"
                          title={isSelf ? "You cannot change your own status" : "Toggle active"}
                        >
                          <Badge tone={u.is_active ? "green" : "slate"}>
                            {u.is_active ? "Active" : "Inactive"}
                          </Badge>
                        </button>
                      </TD>
                      <TD className="text-slate-500">
                        {formatDate(u.created_at)}
                      </TD>
                      <TD className="text-right">
                        <Button
                          size="sm"
                          variant="ghost"
                          disabled={isSelf || busy}
                          onClick={() => setToDelete(u)}
                          aria-label="Delete user"
                        >
                          <Trash2 className="h-4 w-4 text-red-500" />
                        </Button>
                      </TD>
                    </TR>
                  );
                })
              )}
            </TBody>
          </Table>
        </Card>
      )}

      <Modal
        open={!!toDelete}
        onClose={() => setToDelete(null)}
        title="Delete user"
        description="This action cannot be undone."
        footer={
          <>
            <Button variant="outline" onClick={() => setToDelete(null)}>
              Cancel
            </Button>
            <Button variant="danger" onClick={onDelete} loading={deleting}>
              Delete user
            </Button>
          </>
        }
      >
        <div className="flex items-start gap-3">
          <span className="flex h-9 w-9 items-center justify-center rounded-lg bg-red-50 text-red-600">
            <ShieldAlert className="h-5 w-5" />
          </span>
          <p className="text-sm text-slate-600">
            Permanently delete <strong>{toDelete?.email}</strong>? They will lose
            access immediately.
          </p>
        </div>
      </Modal>
    </div>
  );
}
