import type { UploadResult } from "./types";

// Defaults to the gateway's localhost port so `npm run dev` works out of the box.
// In docker-compose, NEXT_PUBLIC_API_URL is baked at build time via the Dockerfile ARG.
export const API_BASE =
  process.env.NEXT_PUBLIC_API_URL || "http://localhost:8080";

export const API_PREFIX = "/api/v1";

const TOKEN_KEY = "g2t_token";

export function getToken(): string | null {
  if (typeof window === "undefined") return null;
  return window.localStorage.getItem(TOKEN_KEY);
}

export function setToken(token: string): void {
  if (typeof window === "undefined") return;
  window.localStorage.setItem(TOKEN_KEY, token);
}

export function clearToken(): void {
  if (typeof window === "undefined") return;
  window.localStorage.removeItem(TOKEN_KEY);
}

export class ApiError extends Error {
  status: number;
  detail: unknown;

  constructor(message: string, status: number, detail: unknown) {
    super(message);
    this.name = "ApiError";
    this.status = status;
    this.detail = detail;
  }
}

function buildUrl(path: string): string {
  // Allow absolute paths; otherwise prefix with /api/v1
  const normalized = path.startsWith("/") ? path : `/${path}`;
  return `${API_BASE}${API_PREFIX}${normalized}`;
}

function authHeaders(extra?: HeadersInit): HeadersInit {
  const headers: Record<string, string> = {};
  const token = getToken();
  if (token) {
    headers["Authorization"] = `Bearer ${token}`;
  }
  if (extra) {
    return { ...headers, ...(extra as Record<string, string>) };
  }
  return headers;
}

function extractDetail(body: unknown): string {
  if (!body) return "";
  if (typeof body === "string") return body;
  if (typeof body === "object") {
    const obj = body as Record<string, unknown>;
    const detail = obj.detail;
    if (typeof detail === "string") return detail;
    if (Array.isArray(detail)) {
      // FastAPI validation error array
      return detail
        .map((d) => {
          if (d && typeof d === "object" && "msg" in d) {
            return String((d as Record<string, unknown>).msg);
          }
          return JSON.stringify(d);
        })
        .join("; ");
    }
    if (detail) return JSON.stringify(detail);
    if (typeof obj.message === "string") return obj.message;
  }
  return "";
}

async function parseResponse<T>(res: Response): Promise<T> {
  const contentType = res.headers.get("content-type") || "";
  let body: unknown = null;
  if (contentType.includes("application/json")) {
    body = await res.json().catch(() => null);
  } else {
    const text = await res.text().catch(() => "");
    body = text || null;
  }

  if (!res.ok) {
    const detail =
      extractDetail(body) || res.statusText || `Request failed (${res.status})`;
    throw new ApiError(detail, res.status, body);
  }

  return body as T;
}

export async function get<T>(path: string): Promise<T> {
  const res = await fetch(buildUrl(path), {
    method: "GET",
    headers: authHeaders(),
    cache: "no-store",
  });
  return parseResponse<T>(res);
}

export async function post<T>(path: string, data?: unknown): Promise<T> {
  const res = await fetch(buildUrl(path), {
    method: "POST",
    headers: authHeaders({ "Content-Type": "application/json" }),
    body: data !== undefined ? JSON.stringify(data) : undefined,
  });
  return parseResponse<T>(res);
}

export async function patch<T>(path: string, data?: unknown): Promise<T> {
  const res = await fetch(buildUrl(path), {
    method: "PATCH",
    headers: authHeaders({ "Content-Type": "application/json" }),
    body: data !== undefined ? JSON.stringify(data) : undefined,
  });
  return parseResponse<T>(res);
}

export async function del<T = void>(path: string): Promise<T> {
  const res = await fetch(buildUrl(path), {
    method: "DELETE",
    headers: authHeaders(),
  });
  return parseResponse<T>(res);
}

export async function uploadFile(file: File): Promise<UploadResult> {
  const form = new FormData();
  form.append("file", file);
  const res = await fetch(buildUrl("/uploads/wsi"), {
    method: "POST",
    headers: authHeaders(), // do NOT set Content-Type; browser sets multipart boundary
    body: form,
  });
  return parseResponse<UploadResult>(res);
}

/**
 * Fetch a protected resource as a Blob (used for authenticated downloads such
 * as report HTML export). Returns the blob and a suggested filename.
 */
export async function fetchBlob(
  path: string
): Promise<{ blob: Blob; filename: string }> {
  const res = await fetch(buildUrl(path), {
    method: "GET",
    headers: authHeaders(),
  });
  if (!res.ok) {
    const text = await res.text().catch(() => "");
    throw new ApiError(
      text || res.statusText || `Download failed (${res.status})`,
      res.status,
      text
    );
  }
  const blob = await res.blob();
  let filename = "download";
  const disposition = res.headers.get("content-disposition");
  if (disposition) {
    const match = /filename\*?=(?:UTF-8'')?"?([^";]+)"?/i.exec(disposition);
    if (match && match[1]) {
      filename = decodeURIComponent(match[1]);
    }
  }
  return { blob, filename };
}

export const api = {
  get,
  post,
  patch,
  del,
  uploadFile,
  fetchBlob,
  getToken,
  setToken,
  clearToken,
};

export default api;
