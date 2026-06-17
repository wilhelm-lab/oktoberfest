# Stage 1: Build Vue 3 frontend
FROM node:22-alpine AS frontend-build
WORKDIR /app/frontend
COPY frontend/package*.json ./
RUN npm ci --silent
COPY frontend/ ./
RUN npm run build

# Stage 2: Python backend + Oktoberfest
FROM python:3.11-slim

# System dependencies
# - mokapot is pip-installable (no binary needed)
# - percolator: NOT included; document as external requirement
# - ThermoRawFileParser: NOT included; document as external requirement (requires Mono + .NET)
#   Users must pre-convert .raw files to .mzML, or install ThermoRawFileParser separately
#   and mount it into the container.
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY backend/pyproject.toml backend/poetry.lock* ./
# Use pip (not poetry) in Docker for simplicity
RUN pip install --no-cache-dir \
    fastapi>=0.111.0 \
    "uvicorn[standard]>=0.29.0" \
    sqlalchemy>=2.0.0 \
    alembic>=1.13.0 \
    "celery[redis]>=5.3.0" \
    redis>=5.0.0 \
    pydantic>=2.7.0 \
    pydantic-settings>=2.3.0 \
    aiofiles>=23.0.0 \
    python-multipart>=0.0.9 \
    oktoberfest>=0.11.0

# Copy backend application code
COPY backend/app/ ./app/
COPY backend/alembic/ ./alembic/
COPY backend/alembic.ini ./

# Copy built frontend into the backend static directory
COPY --from=frontend-build /app/frontend/dist /app/frontend/dist

# Expose API port
EXPOSE 8000

# Default: run uvicorn. Worker uses the same image with a different command.
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
