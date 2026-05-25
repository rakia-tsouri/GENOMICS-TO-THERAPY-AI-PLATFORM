"""WSI / tissue-image upload endpoint.

Saves the file to the shared volume (mounted into both the gateway and the
histopathology service) and returns a path that can be passed as
`wsi_image_path` when creating a job.
"""
import os
import uuid

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile

from ..config import settings
from ..deps import get_current_user
from ..models import User

router = APIRouter(prefix="/uploads", tags=["uploads"])

_ALLOWED_EXT = {".svs", ".tif", ".tiff", ".png", ".jpg", ".jpeg"}


@router.post("/wsi")
async def upload_wsi(file: UploadFile = File(...), _: User = Depends(get_current_user)):
    ext = os.path.splitext(file.filename or "")[1].lower()
    if ext not in _ALLOWED_EXT:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported file type '{ext}'. Allowed: {sorted(_ALLOWED_EXT)}",
        )

    dest_dir = settings.upload_dir
    os.makedirs(dest_dir, exist_ok=True)
    safe_name = f"{uuid.uuid4().hex}{ext}"
    dest_path = os.path.join(dest_dir, safe_name)

    size = 0
    limit = settings.max_upload_mb * 1024 * 1024
    with open(dest_path, "wb") as out:
        while chunk := await file.read(1024 * 1024):
            size += len(chunk)
            if size > limit:
                out.close()
                os.remove(dest_path)
                raise HTTPException(status_code=413, detail="File too large")
            out.write(chunk)

    return {"wsi_image_path": dest_path, "filename": file.filename, "size_bytes": size}
