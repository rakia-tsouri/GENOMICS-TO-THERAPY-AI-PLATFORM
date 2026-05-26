import { Check, Loader2 } from "lucide-react";
import { cn } from "@/lib/utils";
import type { JobStatus } from "@/lib/types";

export interface Step {
  key: string;
  label: string;
}

export const PIPELINE_STEPS: Step[] = [
  { key: "validating", label: "Validating" },
  { key: "analyzing", label: "Analyzing" },
  { key: "predicting_drugs", label: "Drug discovery" },
  { key: "histopathology", label: "Histopathology" },
  { key: "fusion", label: "Cross-modal fusion" },
  { key: "done", label: "Done" },
];

function stepIndex(currentStep: string): number {
  const idx = PIPELINE_STEPS.findIndex((s) => s.key === currentStep);
  return idx; // -1 if not found (e.g. "pending")
}

export default function Stepper({
  currentStep,
  status,
}: {
  currentStep: string;
  status: JobStatus;
}) {
  const activeIdx = stepIndex(currentStep);
  const isFailed = status === "failed";
  const isDone = status === "completed";

  return (
    <ol className="flex flex-col gap-2 sm:flex-row sm:items-start sm:gap-0">
      {PIPELINE_STEPS.map((step, i) => {
        // A step is complete if it's before the active step, or job is done.
        const complete = isDone || (activeIdx > -1 && i < activeIdx);
        const isCurrent = !isDone && !isFailed && i === activeIdx;
        const failedHere = isFailed && i === Math.max(activeIdx, 0);

        let circle: React.ReactNode;
        let circleCls = "border-slate-300 bg-white text-slate-400";
        if (complete) {
          circleCls = "border-green-500 bg-green-500 text-white";
          circle = <Check className="h-4 w-4" />;
        } else if (failedHere) {
          circleCls = "border-red-500 bg-red-500 text-white";
          circle = <span className="text-xs font-bold">!</span>;
        } else if (isCurrent) {
          circleCls = "border-brand-500 bg-brand-50 text-brand-600";
          circle = <Loader2 className="h-4 w-4 animate-spin" />;
        } else {
          circle = <span className="text-xs font-semibold">{i + 1}</span>;
        }

        const lineDone = isDone || (activeIdx > -1 && i < activeIdx);

        return (
          <li key={step.key} className="flex items-center sm:flex-1">
            <div className="flex items-center gap-2 sm:flex-col sm:items-center sm:gap-1.5">
              <span
                className={cn(
                  "flex h-8 w-8 shrink-0 items-center justify-center rounded-full border-2 transition-colors",
                  circleCls
                )}
              >
                {circle}
              </span>
              <span
                className={cn(
                  "text-xs font-medium sm:text-center",
                  isCurrent
                    ? "text-brand-700"
                    : complete
                      ? "text-slate-700"
                      : failedHere
                        ? "text-red-600"
                        : "text-slate-400"
                )}
              >
                {step.label}
              </span>
            </div>
            {i < PIPELINE_STEPS.length - 1 && (
              <div
                className={cn(
                  "mx-2 hidden h-0.5 flex-1 rounded sm:block",
                  lineDone ? "bg-green-500" : "bg-slate-200"
                )}
              />
            )}
          </li>
        );
      })}
    </ol>
  );
}
