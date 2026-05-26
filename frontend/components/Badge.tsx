import { type ReactNode } from "react";
import { cn } from "@/lib/utils";

type Tone =
  | "neutral"
  | "brand"
  | "green"
  | "red"
  | "amber"
  | "blue"
  | "slate";

const toneClasses: Record<Tone, string> = {
  neutral: "bg-slate-100 text-slate-700 border-slate-200",
  brand: "bg-brand-50 text-brand-700 border-brand-200",
  green: "bg-green-50 text-green-700 border-green-200",
  red: "bg-red-50 text-red-700 border-red-200",
  amber: "bg-amber-50 text-amber-700 border-amber-200",
  blue: "bg-blue-50 text-blue-700 border-blue-200",
  slate: "bg-slate-100 text-slate-600 border-slate-200",
};

export default function Badge({
  children,
  tone = "neutral",
  className,
}: {
  children: ReactNode;
  tone?: Tone;
  className?: string;
}) {
  return (
    <span
      className={cn(
        "inline-flex items-center gap-1 rounded-full border px-2 py-0.5 text-xs font-medium",
        toneClasses[tone],
        className
      )}
    >
      {children}
    </span>
  );
}

/** Tone helper for toxicity risk strings. */
export function toxicityTone(risk: string | undefined | null): Tone {
  const r = (risk || "").toLowerCase();
  if (r.includes("low") || r.includes("none")) return "green";
  if (r.includes("med") || r.includes("moderate")) return "amber";
  if (r.includes("high")) return "red";
  return "slate";
}

/** Tone helper for fusion agreement labels. */
export function agreementTone(agreement: string | undefined | null): Tone {
  const a = (agreement || "").toLowerCase();
  if (a.includes("agree") || a.includes("concordant") || a.includes("match"))
    return "green";
  if (a.includes("partial") || a.includes("weak")) return "amber";
  if (a.includes("disagree") || a.includes("discordant") || a.includes("conflict"))
    return "red";
  return "slate";
}
