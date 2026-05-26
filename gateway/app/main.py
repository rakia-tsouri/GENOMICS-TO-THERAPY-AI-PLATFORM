"""API Gateway for the Genomics-to-Therapy AI Platform.

Responsibilities:
  * AuthN/AuthZ (JWT) and user/project management
  * Job lifecycle + orchestration of the 4 microservices
  * Cross-modal fusion and report generation

This is the single entry point the frontend talks to.
"""
import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .database import Base, engine, wait_for_db
from .routers import auth, jobs, projects, reports, uploads, users
from .seed import seed_demo

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger("gateway")

app = FastAPI(title=settings.app_name, version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origin_list,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def on_startup() -> None:
    wait_for_db()  # tolerate Postgres not being ready/resolvable yet
    Base.metadata.create_all(bind=engine)
    seed_demo()
    logger.info("Gateway ready. DB=%s", settings.database_url.split("@")[-1])


for r in (auth.router, users.router, projects.router, jobs.router, reports.router, uploads.router):
    app.include_router(r, prefix=settings.api_prefix)


@app.get("/")
def root():
    return {"service": settings.app_name, "docs": "/docs", "api": settings.api_prefix}


@app.get("/health")
def health():
    return {"status": "healthy"}
