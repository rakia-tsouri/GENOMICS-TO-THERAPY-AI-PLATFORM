"use client";

import Link from "next/link";
import { useRouter } from "next/navigation";
import { useState } from "react";
import { LogOut, Plus, ChevronDown, UserCircle2 } from "lucide-react";
import { useAuth } from "@/lib/auth";
import Button from "./Button";
import Badge from "./Badge";

export default function Topbar() {
  const { user, logout, isAdmin } = useAuth();
  const router = useRouter();
  const [open, setOpen] = useState(false);

  function handleLogout() {
    logout();
    router.replace("/login");
  }

  return (
    <header className="flex h-16 items-center justify-between border-b border-slate-200 bg-white px-4 md:px-6">
      <div className="md:hidden">
        <span className="text-sm font-semibold text-slate-900">
          Genomics→Therapy
        </span>
      </div>

      <div className="hidden md:block" />

      <div className="flex items-center gap-3">
        <Link href="/jobs/new">
          <Button size="sm">
            <Plus className="h-4 w-4" />
            New analysis
          </Button>
        </Link>

        <div className="relative">
          <button
            onClick={() => setOpen((o) => !o)}
            onBlur={() => setTimeout(() => setOpen(false), 150)}
            className="flex items-center gap-2 rounded-lg px-2 py-1.5 text-sm text-slate-700 hover:bg-slate-100"
          >
            <span className="flex h-8 w-8 items-center justify-center rounded-full bg-brand-100 text-brand-700">
              <UserCircle2 className="h-5 w-5" />
            </span>
            <span className="hidden max-w-[160px] truncate sm:inline">
              {user?.email ?? "Account"}
            </span>
            <ChevronDown className="h-4 w-4 text-slate-400" />
          </button>

          {open && (
            <div className="absolute right-0 top-full z-50 mt-1 w-60 rounded-lg border border-slate-200 bg-white py-1 shadow-lg">
              <div className="border-b border-slate-100 px-4 py-3">
                <p className="truncate text-sm font-medium text-slate-900">
                  {user?.full_name || "—"}
                </p>
                <p className="truncate text-xs text-slate-500">{user?.email}</p>
                <div className="mt-2">
                  <Badge tone={isAdmin ? "brand" : "slate"}>
                    {user?.role ?? "researcher"}
                  </Badge>
                </div>
              </div>
              <Link
                href="/settings"
                className="block px-4 py-2 text-sm text-slate-700 hover:bg-slate-50"
              >
                Settings
              </Link>
              <button
                onClick={handleLogout}
                className="flex w-full items-center gap-2 px-4 py-2 text-left text-sm text-red-600 hover:bg-red-50"
              >
                <LogOut className="h-4 w-4" />
                Log out
              </button>
            </div>
          )}
        </div>
      </div>
    </header>
  );
}
