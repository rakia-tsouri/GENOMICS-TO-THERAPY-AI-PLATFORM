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
from .database import Base, SessionLocal, engine
from .models import User
from .security import hash_password
from .routers import auth, jobs, projects, reports, uploads, users

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


def _seed_admin() -> None:
    db = SessionLocal()
    try:
        if db.query(User).count() == 0:
            db.add(
                User(
                    email=settings.admin_email,
                    full_name="Platform Admin",
                    hashed_password=hash_password(settings.admin_password),
                    role="admin",
                )
            )
            db.commit()
            logger.info("Seeded default admin user: %s", settings.admin_email)
    finally:
        db.close()


@app.on_event("startup")
def on_startup() -> None:
    Base.metadata.create_all(bind=engine)
    _seed_admin()
    logger.info("Gateway ready. DB=%s", settings.database_url.split("@")[-1])


for r in (auth.router, users.router, projects.router, jobs.router, reports.router, uploads.router):
    app.include_router(r, prefix=settings.api_prefix)


@app.get("/")
def root():
    return {"service": settings.app_name, "docs": "/docs", "api": settings.api_prefix}


@app.get("/health")
def health():
    return {"status": "healthy"}
