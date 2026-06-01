"use client";

import { useEffect, useRef, useState } from "react";

/**
 * Render a SMILES string as a 2D molecular structure using smiles-drawer.
 * The library is dynamically imported so it only ships in the bundle
 * where it is used.
 */
export default function MoleculeImage({
  smiles,
  width = 180,
  height = 120,
}: {
  smiles: string;
  width?: number;
  height?: number;
}) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [failed, setFailed] = useState(false);

  useEffect(() => {
    let cancelled = false;

    type DrawerInstance = {
      draw: (parsed: unknown, target: HTMLCanvasElement, theme: string) => void;
    };
    type SDApi = {
      Drawer: new (opts: Record<string, unknown>) => DrawerInstance;
      parse: (
        s: string,
        ok: (tree: unknown) => void,
        err?: (e: unknown) => void
      ) => void;
    };

    async function draw() {
      if (!canvasRef.current || !smiles) return;
      try {
        // smiles-drawer is a UMD module; resolve the namespace at runtime.
        const mod = (await import("smiles-drawer")) as unknown as
          | (SDApi & { default?: SDApi })
          | { default: SDApi };
        const sd: SDApi =
          (mod as { default?: SDApi }).default ?? (mod as SDApi);
        if (cancelled || !canvasRef.current) return;
        const drawer = new sd.Drawer({
          width,
          height,
          padding: 8,
          bondThickness: 1.2,
        });
        sd.parse(
          smiles,
          (tree) => {
            if (!cancelled && canvasRef.current) {
              drawer.draw(tree, canvasRef.current, "light");
            }
          },
          () => setFailed(true)
        );
      } catch {
        setFailed(true);
      }
    }
    void draw();
    return () => {
      cancelled = true;
    };
  }, [smiles, width, height]);

  if (failed) {
    return (
      <div
        className="flex items-center justify-center rounded border border-slate-200 bg-slate-50 text-[10px] text-slate-400"
        style={{ width, height }}
      >
        SMILES invalide
      </div>
    );
  }

  return (
    <canvas
      ref={canvasRef}
      width={width}
      height={height}
      className="rounded border border-slate-200 bg-white"
    />
  );
}
