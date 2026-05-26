import { cn } from "@/lib/utils";
import type { JobStatus } from "@/lib/types";

const config: Record<
  JobStatus,
  { label: string; className: string; dot: string; pulse?: boolean }
> = {
  pending: {
    label: "Pending",
    className: "bg-slate-100 text-slate-600 border-slate-200",
    dot: "bg-slate-400",
  },
  running: {
    label: "Running",
    className: "bg-blue-50 text-blue-700 border-blue-200",
    dot: "bg-blue-500",
    pulse: true,
  },
  completed: {
    label: "Completed",
    className: "bg-green-50 text-green-700 border-green-200",
    dot: "bg-green-500",
  },
  failed: {
    label: "Failed",
    className: "bg-red-50 text-red-700 border-red-200",
    dot: "bg-red-500",
  },
};

export default function StatusBadge({
  status,
  className,
}: {
  status: JobStatus;
  className?: string;
}) {
  const c = config[status] ?? config.pending;
  return (
    <span
      className={cn(
        "inline-flex items-center gap-1.5 rounded-full border px-2.5 py-0.5 text-xs font-medium",
        c.className,
        className
      )}
    >
      <span
        className={cn(
          "h-1.5 w-1.5 rounded-full",
          c.dot,
          c.pulse && "animate-pulse"
        )}
      />
      {c.label}
    </span>
  );
}
