from __future__ import annotations

from fastapi import APIRouter
from celery.app.control import Inspect

from app.config import settings

router = APIRouter(tags=["health"])


@router.get("/health")
def health_check():
    worker_status = "unknown"
    try:
        from app.worker.celery_app import celery_app

        insp = celery_app.control.inspect(timeout=1.0)
        active = insp.ping()
        worker_status = "up" if active else "down"
    except Exception:
        worker_status = "down"

    return {"status": "ok", "mode": settings.app_mode, "worker": worker_status}
