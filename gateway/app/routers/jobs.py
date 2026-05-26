from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from ..database import get_db
from ..deps import get_current_user
from ..models import Job, User
from ..orchestration.pipeline import run_pipeline
from ..schemas.job import JobCreate, JobRead, JobSummary

router = APIRouter(prefix="/jobs", tags=["jobs"])


def _get_owned(db: Session, job_id: int, user: User) -> Job:
    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.owner_id != user.id and user.role != "admin":
        raise HTTPException(status_code=403, detail="Not your job")
    return job


@router.post("", response_model=JobRead, status_code=201)
def create_job(
    body: JobCreate,
    background: BackgroundTasks,
    user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    """Create a pipeline run and kick it off in the background.

    Returns immediately with status='pending'; poll GET /jobs/{id} for progress.
    """
    job = Job(
        name=body.name,
        owner_id=user.id,
        project_id=body.project_id,
        dna_sequence=body.dna_sequence,
        gene_id=body.gene_id,
        wsi_image_path=body.wsi_image_path,
        target_mutations=body.target_mutations,
        status="pending",
        current_step="queued",
    )
    db.add(job)
    db.commit()
    db.refresh(job)

    # run_pipeline is async; BackgroundTasks awaits coroutine functions.
    background.add_task(run_pipeline, job.id)
    return job


@router.get("", response_model=list[JobSummary])
def list_jobs(
    project_id: int | None = Query(default=None),
    user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    q = db.query(Job).filter(Job.owner_id == user.id)
    if project_id is not None:
        q = q.filter(Job.project_id == project_id)
    return q.order_by(Job.created_at.desc()).all()


@router.get("/{job_id}", response_model=JobRead)
def get_job(job_id: int, user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    return _get_owned(db, job_id, user)


@router.delete("/{job_id}", status_code=204)
def delete_job(job_id: int, user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    job = _get_owned(db, job_id, user)
    db.delete(job)
    db.commit()
