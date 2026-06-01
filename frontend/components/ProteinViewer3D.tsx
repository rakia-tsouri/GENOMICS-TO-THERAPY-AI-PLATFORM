"use client";

import { useEffect, useRef, useState } from "react";

declare global {
  interface Window {
    $3Dmol?: {
      createViewer: (el: HTMLElement, opts?: Record<string, unknown>) => Viewer3D;
      SurfaceType: { VDW: number };
    };
  }
}

type Viewer3D = {
  addModel: (data: string, format: string) => void;
  setStyle: (sel: Record<string, unknown>, style: Record<string, unknown>) => void;
  addSurface: (type: number, opts: Record<string, unknown>) => void;
  zoomTo: () => void;
  render: () => void;
  spin: (yes: boolean) => void;
  clear: () => void;
};

export default function ProteinViewer3D({
  pdbPath,
  height = "h-72",
}: {
  pdbPath: string;
  height?: string;
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [error, setError] = useState<string | null>(null);
  const [ready, setReady] = useState(false);

  useEffect(() => {
    let cancelled = false;
    let viewer: Viewer3D | null = null;

    async function init() {
      const $3Dmol = window.$3Dmol;
      if (!$3Dmol || !containerRef.current) return;
      try {
        const text = await fetch(pdbPath).then((r) => {
          if (!r.ok) throw new Error(`PDB fetch failed (${r.status})`);
          return r.text();
        });
        if (cancelled || !containerRef.current) return;
        containerRef.current.innerHTML = "";
        viewer = $3Dmol.createViewer(containerRef.current, {
          backgroundColor: "0x0f172a",
        });
        viewer.addModel(text, "pdb");
        viewer.setStyle({}, { cartoon: { colorscheme: "spectrum" } });
        // Render ligands as stick if any
        viewer.setStyle({ hetflag: true }, {
          stick: { colorscheme: "Jmol", radius: 0.18 },
        });
        viewer.zoomTo();
        viewer.render();
        viewer.spin(true);
        setReady(true);
      } catch (e) {
        setError(e instanceof Error ? e.message : "Failed to load structure");
      }
    }

    // wait for 3Dmol to be on window, then init
    const iv = setInterval(() => {
      if (window.$3Dmol) {
        clearInterval(iv);
        void init();
      }
    }, 100);

    return () => {
      cancelled = true;
      clearInterval(iv);
      try {
        viewer?.spin(false);
        viewer?.clear();
      } catch {}
    };
  }, [pdbPath]);

  return (
    <div className="relative w-full overflow-hidden rounded-lg border border-slate-200 bg-slate-900">
      <div ref={containerRef} className={`${height} w-full`} />
      {!ready && !error && (
        <div className="absolute inset-0 flex items-center justify-center text-xs text-slate-300">
          Chargement du modèle 3D…
        </div>
      )}
      {error && (
        <div className="absolute inset-0 flex items-center justify-center px-4 text-center text-xs text-amber-300">
          {error}
        </div>
      )}
    </div>
  );
}
