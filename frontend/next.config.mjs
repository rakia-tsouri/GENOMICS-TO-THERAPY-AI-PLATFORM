/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  // smiles-drawer 2.3.x ships TypeScript sources in node_modules; let Next
  // transpile it so 'next build' doesn't fail on the raw .ts files.
  transpilePackages: ["smiles-drawer"],
};

export default nextConfig;
