import { type ReactNode } from "react";
import { cn } from "@/lib/utils";

export interface StatChipProps {
  label: string;
  value: ReactNode;
  icon?: ReactNode;
  hint?: string;
  className?: string;
}

export default function StatChip({
  label,
  value,
  icon,
  hint,
  className,
}: StatChipProps) {
  return (
    <div
      className={cn(
        "rounded-xl border border-slate-200 bg-white p-5 shadow-sm",
        className
      )}
    >
      <div className="flex items-center justify-between">
        <p className="text-xs font-medium uppercase tracking-wide text-slate-500">
          {label}
        </p>
        {icon && (
          <span className="flex h-8 w-8 items-center justify-center rounded-lg bg-brand-50 text-brand-600">
            {icon}
          </span>
        )}
      </div>
      <p className="mt-3 text-2xl font-semibold text-slate-900">{value}</p>
      {hint && <p className="mt-1 text-xs text-slate-500">{hint}</p>}
    </div>
  );
}
