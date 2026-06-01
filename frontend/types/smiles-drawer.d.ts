// smiles-drawer ships no .d.ts; the typed shape we need lives in
// components/MoleculeImage.tsx — declare the module as opaque here so TS
// stops complaining and lets the runtime cast in MoleculeImage do its job.
declare module "smiles-drawer";
