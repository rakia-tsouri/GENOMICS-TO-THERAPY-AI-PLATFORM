from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import HTMLResponse
from sqlalchemy.orm import Session

from ..database import get_db
from ..deps import get_current_user
from ..models import Job, Report, User
from ..schemas.report import ReportCreate, ReportRead
from ..services.report_builder import build_summary, render_html

router = APIRouter(prefix="/reports", tags=["reports"])


def _get_owned(db: Session, report_id: int, user: User) -> Report:
    report = db.query(Report).filter(Report.id == report_id).first()
    if not report:
        raise HTTPException(status_code=404, detail="Report not found")
    if report.owner_id != user.id and user.role != "admin":
        raise HTTPException(status_code=403, detail="Not your report")
    return report


@router.post("", response_model=ReportRead, status_code=201)
def create_report(body: ReportCreate, user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    job = db.query(Job).filter(Job.id == body.job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.owner_id != user.id and user.role != "admin":
        raise HTTPException(status_code=403, detail="Not your job")
    if job.status != "completed":
        raise HTTPException(status_code=400, detail="Job is not completed yet")

    summary = build_summary(job)
    report = Report(
        job_id=job.id,
        owner_id=user.id,
        title=body.title or f"Report — {job.name}",
        summary=summary,
    )
    db.add(report)
    db.commit()
    db.refresh(report)
    return report


@router.get("", response_model=list[ReportRead])
def list_reports(user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    return (
        db.query(Report)
        .filter(Report.owner_id == user.id)
        .order_by(Report.created_at.desc())
        .all()
    )


@router.get("/{report_id}", response_model=ReportRead)
def get_report(report_id: int, user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    return _get_owned(db, report_id, user)


@router.get("/{report_id}/export", response_class=HTMLResponse)
def export_report(report_id: int, user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Render the report as a standalone HTML document (printable to PDF)."""
    report = _get_owned(db, report_id, user)
    return HTMLResponse(render_html(report.summary or {}))


@router.delete("/{report_id}", status_code=204)
def delete_report(report_id: int, user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    report = _get_owned(db, report_id, user)
    db.delete(report)
    db.commit()
