# Genomics→Therapy AI Platform — Frontend

A Next.js (App Router) + TypeScript + Tailwind CSS dashboard for the
**Genomics-to-Therapy AI Platform**. It is the UI for an existing FastAPI
gateway and covers authentication, projects, multi-stage analyses (jobs),
results, reports, and admin user management.

## Stack

- **Next.js 14.2** (App Router) + **React 18**
- **TypeScript 5**
- **Tailwind CSS 3.4** (indigo/violet accent `#6c5ce7`)
- `lucide-react` icons, `clsx` for class merging
- JWT bearer auth, token stored in `localStorage`

## Getting started

> The repo ships **without** `node_modules`. Install dependencies first.

```bash
cd frontend
npm install
cp .env.example .env.local   # set NEXT_PUBLIC_API_URL if not localhost:8080
npm run dev
```

Then open http://localhost:3000.

## Environment

| Variable              | Default                 | Description                       |
| --------------------- | ----------------------- | --------------------------------- |
| `NEXT_PUBLIC_API_URL` | `http://localhost:8080` | Base URL of the FastAPI gateway.  |

All API calls target `${NEXT_PUBLIC_API_URL}/api/v1`.

## Scripts

| Command         | Description                |
| --------------- | -------------------------- |
| `npm run dev`   | Start the dev server       |
| `npm run build` | Production build           |
| `npm start`     | Serve the production build |
| `npm run lint`  | Run ESLint                 |

## Project structure

```
app/
  layout.tsx              Root layout (AuthProvider)
  page.tsx                Redirects to /dashboard
  globals.css
  login/ , register/      Public auth pages
  (app)/                  Authenticated area (AppShell + auth guard)
    layout.tsx
    dashboard/            Stat cards + recent analyses
    jobs/                 List, /new creation form, /[id] results
    projects/            List + create, /[id] detail
    reports/             List, /[id] detail (+ authenticated HTML export)
    admin/users/         Admin-only user management
    settings/            Profile & password
components/               Hand-rolled UI primitives + layout
lib/
  api.ts                 fetch wrapper (get/post/patch/del/uploadFile/fetchBlob)
  auth.tsx               AuthProvider / useAuth context
  types.ts               Shared API types
  utils.ts               Formatting + cn() helpers
```

## Key behaviors

- **Auth guard**: `AppShell` shows a spinner while loading, then redirects to
  `/login` if there is no authenticated user.
- **Job results polling**: the job detail page polls `GET /jobs/{id}` every 3s
  while the job is `pending`/`running`, drives a pipeline stepper from
  `current_step`, and renders stage result cards on completion.
- **Authenticated report export**: report HTML export requires the bearer token,
  so it is fetched as a blob and opened in a new tab (`fetchBlob`).
- **Admin**: the Users nav link and page are only shown/usable to `admin` roles.

## Docker

```bash
docker build -t g2t-frontend .
docker run -p 3000:3000 -e NEXT_PUBLIC_API_URL=http://host.docker.internal:8080 g2t-frontend
```

> For research use only. Predictions are not clinical advice.
