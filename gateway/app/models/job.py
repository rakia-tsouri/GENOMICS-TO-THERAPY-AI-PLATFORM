from datetime import datetime, timezone

from sqlalchemy import DateTime, ForeignKey, String, Text, JSON
from sqlalchemy.orm import Mapped, mapped_column, relationship

from ..database import Base


def _now() -> datetime:
    return datetime.now(timezone.utc)


class Job(Base):
    """A single end-to-end pipeline run (DNA and/or WSI -> drug candidates)."""

    __tablename__ = "jobs"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(255), default="Untitled analysis")
    owner_id: Mapped[int] = mapped_column(ForeignKey("users.id"), index=True)
    project_id: Mapped[int | None] = mapped_column(ForeignKey("projects.id"), nullable=True, index=True)

    # Lifecycle: pending -> running -> completed | failed
    status: Mapped[str] = mapped_column(String(32), default="pending", index=True)
    current_step: Mapped[str] = mapped_column(String(64), default="queued")
    error: Mapped[str] = mapped_column(Text, default="")

    # Inputs
    dna_sequence: Mapped[str | None] = mapped_column(Text, nullable=True)
    gene_id: Mapped[str | None] = mapped_column(String(128), nullable=True)
    wsi_image_path: Mapped[str | None] = mapped_column(String(512), nullable=True)
    target_mutations: Mapped[list | None] = mapped_column(JSON, nullable=True)

    # Per-stage outputs (raw JSON from each microservice) + fusion
    validation_result: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    analysis_result: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    drug_result: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    histopathology_result: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    fusion_result: Mapped[dict | None] = mapped_column(JSON, nullable=True)

    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=_now)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)

    owner = relationship("User", back_populates="jobs")
    project = relationship("Project", back_populates="jobs")
    reports = relationship("Report", back_populates="job", cascade="all, delete-orphan")
