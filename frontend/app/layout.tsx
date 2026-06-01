import type { Metadata } from "next";
import Script from "next/script";
import "./globals.css";
import { AuthProvider } from "@/lib/auth";

export const metadata: Metadata = {
  title: "Genomics→Therapy AI Platform",
  description:
    "End-to-end genomics, structure, drug discovery and histopathology AI pipeline.",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body>
        <AuthProvider>{children}</AuthProvider>
        {/* 3Dmol.js for the protein structure viewer (small, loads from CDN) */}
        <Script src="https://3dmol.org/build/3Dmol-min.js" strategy="afterInteractive" />
      </body>
    </html>
  );
}
