from __future__ import annotations

from celery import Celery

from app.config import settings

celery_app = Celery(
    "oktoberfest_webui",
    broker=settings.redis_url,
    backend=settings.redis_url,
    include=["app.worker.tasks"],
)

celery_app.conf.update(
    task_acks_late=True,
    worker_prefetch_multiplier=1,
    task_serializer="json",
    result_serializer="json",
    accept_content=["json"],
    timezone="UTC",
    enable_utc=True,
    # Generous time limits — proteomics jobs can be very long
    task_soft_time_limit=86400,  # 24h soft
    task_time_limit=90000,  # 25h hard
    worker_concurrency=settings.celery_concurrency,
    # task_routes={"run_oktoberfest_job": {"queue": "heavy"}},
)

from celery.signals import worker_process_init
from app.db import engine

@worker_process_init.connect
def dispose_sqlalchemy_engine(**kwargs):
    """Dispose of the connection pool after Celery forks a worker process.
    This prevents 'sqlite3.OperationalError: unable to open database file'
    caused by sharing SQLite connections across fork boundaries."""
    engine.dispose()
