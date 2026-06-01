"""Idempotent demo seeding.

Creates two demo accounts (whose credentials are shown on the login page),
four realistic projects, five completed analyses with curated results, and a
few pre-generated reports. The whole workspace looks populated on first run.
"""
import logging
from datetime import datetime, timedelta, timezone

from sqlalchemy.orm import Session

from .config import settings
from .database import SessionLocal
from .demo_data import JOBS as DEMO_JOBS, PROJECTS as DEMO_PROJECTS, REPORT_FOR_JOBS
from .models import Job, Project, Report, User
from .security import hash_password
from .services.report_builder import build_summary

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


def _ensure_projects(db: Session, owner: User) -> dict[str, Project]:
    by_name: dict[str, Project] = {}
    for spec in DEMO_PROJECTS:
        existing = (
            db.query(Project)
            .filter(Project.owner_id == owner.id, Project.name == spec["name"])
            .first()
        )
        if existing:
            by_name[spec["name"]] = existing
            continue
        proj = Project(owner_id=owner.id, **spec)
        db.add(proj)
        db.commit()
        db.refresh(proj)
        by_name[spec["name"]] = proj
        logger.info("Seeded project: %s", proj.name)
    return by_name


def _ensure_jobs(db: Session, owner: User, projects: dict[str, Project]) -> list[Job]:
    inserted: list[Job] = []
    now = datetime.now(timezone.utc)
    for spec in DEMO_JOBS:
        if (
            db.query(Job)
            .filter(Job.owner_id == owner.id, Job.name == spec["name"])
            .first()
        ):
            continue
        days_ago = spec.get("days_ago", 1)
        created = now - timedelta(days=days_ago)
        job = Job(
            name=spec["name"],
            owner_id=owner.id,
            project_id=projects[spec["project"]].id,
            status="completed",
            current_step="done",
            error="",
            gene_id=spec.get("gene_id"),
            dna_sequence=spec.get("dna_sequence"),
            wsi_image_path=spec.get("wsi_image_path"),
            target_mutations=spec.get("target_mutations"),
            validation_result=spec.get("validation_result"),
            analysis_result=spec.get("analysis_result"),
            drug_result=spec.get("drug_result"),
            histopathology_result=spec.get("histopathology_result"),
            fusion_result=spec.get("fusion_result"),
            created_at=created,
            completed_at=created + timedelta(minutes=4),
        )
        db.add(job)
        db.commit()
        db.refresh(job)
        inserted.append(job)
        logger.info("Seeded job: %s", job.name)
    return inserted


def _ensure_reports(db: Session, owner: User) -> None:
    for job_name in REPORT_FOR_JOBS:
        job = (
            db.query(Job)
            .filter(Job.owner_id == owner.id, Job.name == job_name)
            .first()
        )
        if not job:
            continue
        title = f"Report — {job.name}"
        if (
            db.query(Report)
            .filter(Report.job_id == job.id, Report.title == title)
            .first()
        ):
            continue
        db.add(Report(
            job_id=job.id,
            owner_id=owner.id,
            title=title,
            summary=build_summary(job),
        ))
        db.commit()
        logger.info("Seeded report for: %s", job.name)


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
        projects = _ensure_projects(db, researcher)
        _ensure_jobs(db, researcher, projects)
        _ensure_reports(db, researcher)
    finally:
        db.close()
