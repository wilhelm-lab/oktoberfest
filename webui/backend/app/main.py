from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from app.config import settings
from app.db import Base, engine
import app.models  # noqa: F401 — registers Job model with Base.metadata before create_all

# Ensure data dir exists
Path(settings.data_dir).mkdir(parents=True, exist_ok=True)

# Create DB tables (Alembic handles migrations; this is a safety net for fresh installs)
Base.metadata.create_all(bind=engine)

app = FastAPI(
    title="Oktoberfest Web UI",
    version="0.1.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc",
    openapi_url="/api/openapi.json",
)

# CORS
origins = [o.strip() for o in settings.cors_origins.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins if origins != ["*"] else ["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# API routers
from app.api.v1 import jobs, files, meta, health  # noqa: E402

app.include_router(health.router, prefix="/api/v1")
app.include_router(meta.router, prefix="/api/v1/meta")
app.include_router(jobs.router, prefix="/api/v1/jobs")
app.include_router(files.router, prefix="/api/v1/jobs")

# Serve built SPA (frontend/dist) — check Docker path (/app/frontend/dist) and
# local-dev path (../../frontend/dist relative to app/main.py).
_here = Path(__file__).parent
_dist_candidates = [
    _here.parent / "frontend" / "dist",          # Docker: WORKDIR=/app, app at /app/app/
    _here.parent.parent / "frontend" / "dist",   # local dev: webui/backend/app/
]
_dist = next((p for p in _dist_candidates if p.exists()), None)

if _dist is not None:
    app.mount("/assets", StaticFiles(directory=str(_dist / "assets")), name="assets")

    from fastapi.responses import FileResponse

    @app.get("/{full_path:path}", include_in_schema=False)
    async def spa_fallback(full_path: str):
        index = _dist / "index.html"
        return FileResponse(str(index))
