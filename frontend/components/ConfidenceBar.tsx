import { cn, toPercentValue } from "@/lib/utils";

export interface ConfidenceBarProps {
  /** A value either 0-1 (fraction) or 0-100 (percent). */
  value: number | null | undefined;
  label?: string;
  showValue?: boolean;
  className?: string;
  /** Override the auto color thresholds. */
  color?: string;
}

function autoColor(pct: number): string {
  if (pct >= 75) return "bg-green-500";
  if (pct >= 50) return "bg-amber-500";
  if (pct >= 25) return "bg-orange-500";
  return "bg-red-500";
}

export default function ConfidenceBar({
  value,
  label,
  showValue = true,
  className,
  color,
}: ConfidenceBarProps) {
  const pct = toPercentValue(value);
  const barColor = color ?? autoColor(pct);

  return (
    <div className={cn("w-full", className)}>
      {(label || showValue) && (
        <div className="mb-1 flex items-center justify-between text-xs text-slate-600">
          {label && <span>{label}</span>}
          {showValue && (
            <span className="font-medium tabular-nums">{pct.toFixed(0)}%</span>
          )}
        </div>
      )}
      <div className="h-2 w-full overflow-hidden rounded-full bg-slate-100">
        <div
          className={cn("h-full rounded-full transition-all", barColor)}
          style={{ width: `${pct}%` }}
        />
      </div>
    </div>
  );
}
