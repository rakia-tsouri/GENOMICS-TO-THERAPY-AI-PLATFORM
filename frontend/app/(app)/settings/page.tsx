"use client";

import { useEffect, useState } from "react";
import { UserCog, KeyRound } from "lucide-react";
import { patch, ApiError } from "@/lib/api";
import type { UpdateMeInput, User } from "@/lib/types";
import { useAuth } from "@/lib/auth";
import PageHeader from "@/components/PageHeader";
import Card, { CardHeader, CardBody } from "@/components/Card";
import Button from "@/components/Button";
import Badge from "@/components/Badge";
import { Field, Input } from "@/components/Input";
import {
  ErrorBanner,
  SuccessBanner,
  PageSpinner,
} from "@/components/Feedback";

export default function SettingsPage() {
  const { user, setUser, loading } = useAuth();

  const [fullName, setFullName] = useState("");
  const [savingProfile, setSavingProfile] = useState(false);
  const [profileError, setProfileError] = useState<string | null>(null);
  const [profileOk, setProfileOk] = useState<string | null>(null);

  const [password, setPassword] = useState("");
  const [confirm, setConfirm] = useState("");
  const [savingPw, setSavingPw] = useState(false);
  const [pwError, setPwError] = useState<string | null>(null);
  const [pwOk, setPwOk] = useState<string | null>(null);

  useEffect(() => {
    if (user) setFullName(user.full_name);
  }, [user]);

  async function onSaveProfile(e: React.FormEvent) {
    e.preventDefault();
    setProfileError(null);
    setProfileOk(null);
    setSavingProfile(true);
    try {
      const payload: UpdateMeInput = { full_name: fullName };
      const updated = await patch<User>("/users/me", payload);
      setUser(updated);
      setProfileOk("Profile updated.");
    } catch (err) {
      setProfileError(
        err instanceof ApiError ? err.message : "Failed to update profile."
      );
    } finally {
      setSavingProfile(false);
    }
  }

  async function onChangePassword(e: React.FormEvent) {
    e.preventDefault();
    setPwError(null);
    setPwOk(null);
    if (password.length < 8) {
      setPwError("Password must be at least 8 characters.");
      return;
    }
    if (password !== confirm) {
      setPwError("Passwords do not match.");
      return;
    }
    setSavingPw(true);
    try {
      const payload: UpdateMeInput = { password };
      await patch<User>("/users/me", payload);
      setPassword("");
      setConfirm("");
      setPwOk("Password changed successfully.");
    } catch (err) {
      setPwError(
        err instanceof ApiError ? err.message : "Failed to change password."
      );
    } finally {
      setSavingPw(false);
    }
  }

  if (loading || !user) return <PageSpinner />;

  return (
    <div>
      <PageHeader
        title="Settings"
        description="Manage your profile and credentials."
      />

      <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
        <Card>
          <CardHeader
            title="Profile"
            description="Your account information."
            icon={<UserCog className="h-4 w-4" />}
          />
          <CardBody>
            <form onSubmit={onSaveProfile} className="space-y-4">
              <SuccessBanner message={profileOk} />
              <ErrorBanner message={profileError} />

              <Field label="Email">
                <Input value={user.email} disabled />
              </Field>

              <Field label="Role">
                <div>
                  <Badge tone={user.role === "admin" ? "brand" : "slate"}>
                    {user.role}
                  </Badge>
                </div>
              </Field>

              <Field label="Full name" htmlFor="full_name" required>
                <Input
                  id="full_name"
                  required
                  value={fullName}
                  onChange={(e) => setFullName(e.target.value)}
                />
              </Field>

              <div className="flex justify-end">
                <Button
                  type="submit"
                  loading={savingProfile}
                  disabled={fullName.trim() === user.full_name || !fullName.trim()}
                >
                  Save changes
                </Button>
              </div>
            </form>
          </CardBody>
        </Card>

        <Card>
          <CardHeader
            title="Change password"
            description="Use a strong, unique password."
            icon={<KeyRound className="h-4 w-4" />}
          />
          <CardBody>
            <form onSubmit={onChangePassword} className="space-y-4">
              <SuccessBanner message={pwOk} />
              <ErrorBanner message={pwError} />

              <Field
                label="New password"
                htmlFor="new_pw"
                hint="At least 8 characters."
                required
              >
                <Input
                  id="new_pw"
                  type="password"
                  autoComplete="new-password"
                  required
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                />
              </Field>

              <Field label="Confirm new password" htmlFor="confirm_pw" required>
                <Input
                  id="confirm_pw"
                  type="password"
                  autoComplete="new-password"
                  required
                  value={confirm}
                  onChange={(e) => setConfirm(e.target.value)}
                />
              </Field>

              <div className="flex justify-end">
                <Button
                  type="submit"
                  loading={savingPw}
                  disabled={!password || !confirm}
                >
                  Update password
                </Button>
              </div>
            </form>
          </CardBody>
        </Card>
      </div>
    </div>
  );
}
