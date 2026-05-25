"""Idempotent demo seeding.

Creates two demo accounts (whose credentials are shown on the login page) and a
sample project, so the platform is usable immediately on a fresh database.
Safe to run on every startup: existing rows are left untouched.
"""
import logging

from sqlalchemy.orm import Session

from .config import settings
from .database import SessionLocal
from .models import Project, User
from .security import hash_password

logger = logging.getLogger("gateway.seed")


def _ensure_user(db: Session, email: str, password: str, full_name: str, role: str) -> User:
    user = db.query(User).filter(User.email == email).first()
    if user:
        return user
    user = User(
        email=email,
        full_name=full_name,
        hashed_password=hash_password(password),
        role=role,
    )
    db.add(user)
    db.commit()
    db.refresh(user)
    logger.info("Seeded %s account: %s", role, email)
    return user


def seed_demo() -> None:
    if not settings.seed_demo:
        return
    db = SessionLocal()
    try:
        _ensure_user(
            db, settings.demo_admin_email, settings.demo_admin_password,
            "Demo Admin", "admin",
        )
        researcher = _ensure_user(
            db, settings.demo_researcher_email, settings.demo_researcher_password,
            "Demo Researcher", "researcher",
        )

        # A sample project so the workspace isn't empty on first login.
        sample_name = "TP53 — Lung Adenocarcinoma (sample)"
        exists = (
            db.query(Project)
            .filter(Project.owner_id == researcher.id, Project.name == sample_name)
            .first()
        )
        if not exists:
            db.add(Project(
                name=sample_name,
                description=(
                    "Sample project. Try an analysis with gene ID NM_000546 (TP53) "
                    "or upload a tissue image for the Track B demo."
                ),
                cancer_type="LUAD",
                owner_id=researcher.id,
            ))
            db.commit()
            logger.info("Seeded sample project for %s", researcher.email)
    finally:
        db.close()
