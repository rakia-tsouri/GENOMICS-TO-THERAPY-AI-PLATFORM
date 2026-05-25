from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel, ConfigDict


class ReportCreate(BaseModel):
    job_id: int
    title: Optional[str] = None


class ReportRead(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    job_id: int
    owner_id: int
    title: str
    summary: Optional[dict[str, Any]]
    created_at: datetime
