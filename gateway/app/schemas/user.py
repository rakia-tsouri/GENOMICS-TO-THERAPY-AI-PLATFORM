from datetime import datetime
from typing import Optional

from pydantic import BaseModel, ConfigDict, EmailStr


class UserRead(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    email: EmailStr
    full_name: str
    role: str
    is_active: bool
    created_at: datetime


class UserUpdate(BaseModel):
    full_name: Optional[str] = None
    password: Optional[str] = None


class AdminUserUpdate(BaseModel):
    """Fields an admin may change on any user."""
    full_name: Optional[str] = None
    role: Optional[str] = None       # "researcher" | "admin"
    is_active: Optional[bool] = None
