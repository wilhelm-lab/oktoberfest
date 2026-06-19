# Oktoberfest Web UI

A browser-based interface for [Oktoberfest](https://github.com/wilhelm-lab/oktoberfest) — open-source spectral library generation and rescoring pipeline based on Prosit.

Submit and monitor **Rescoring**, **Collision Energy Calibration**, and **Spectral Library Generation** jobs through a guided form UI, without writing JSON configs by hand.

---

## Quickstart (Docker Compose)

```bash
# 1. Clone and enter the webui directory
cd oktoberfest/webui

# 2. Copy environment file
cp .env.example .env

# 3. Build and start
docker compose up --build

# 4. Open http://localhost:8000
```

The `api` and `worker` services share a persistent `okt-data` Docker volume. Your jobs survive container restarts.

---

## Local development (no Docker)

Prerequisites: Python 3.10+, Node 18+, Redis, `oktoberfest` installed, and a local Python environment of your choice.

Create and activate a Python environment first if you do not already have one. Then install poetry and the backend dependencies:

```bash
pip install poetry
```

```bash
cd backend
poetry install
```

From there, the backend can be run with Poetry-managed commands:

```bash
# 1. Database
poetry run alembic upgrade head

# 2. Start Redis (separate terminal or Docker)
docker run -p 6379:6379 redis:7-alpine

# 3. API (separate terminal)
poetry run uvicorn app.main:app --reload --port 8000

# 4. Worker (separate terminal)
poetry run celery -A app.worker.celery_app worker --loglevel=info --concurrency=1

# 5. Frontend dev server (separate terminal)
cd ../frontend && npm install && npm run dev
# → http://localhost:5173  (proxies /api → :8000)
```

---

## Architecture

```
Browser (Vue 3 SPA)
  │  POST /api/v1/jobs  (create + upload + submit)
  │  GET  /api/v1/jobs/{id}  (poll status every ~3s)
  ▼
FastAPI backend (port 8000)
  │  enqueue task  (Celery + Redis)
  ▼
Celery worker
  │  python -m oktoberfest --config_path config.json
  │  zip output/ → results.zip
  ▼
SQLite (./data/app.db) + DATA_DIR (./data/<job_id>/)
```

Each job gets its own directory:

```
data/<job_id>/
  inputs/        ← uploaded files
  output/        ← Oktoberfest writes here
  config.json    ← resolved config (exact copy in results.zip)
  captured.log   ← stdout/stderr from oktoberfest process
  results.zip    ← downloadable archive
```

---

## Job types

| Type                 | Required inputs          | Outputs                                |
| -------------------- | ------------------------ | -------------------------------------- |
| **Rescoring**        | Search results + spectra | Rescored PSM/peptide tables, FDR plots |
| **CE Calibration**   | Search results + spectra | Optimal NCE + diagnostics              |
| **Spectral Library** | FASTA or peptides CSV    | `.msp` / Spectronaut / DLib library    |

---

## Configuration

All settings are environment variables (see `.env.example`):

| Variable             | Default                    | Description                         |
| -------------------- | -------------------------- | ----------------------------------- |
| `DATA_DIR`           | `./data`                   | Job storage root                    |
| `DATABASE_URL`       | `sqlite:///./data/app.db`  | SQLite (local) or Postgres URL      |
| `REDIS_URL`          | `redis://localhost:6379/0` | Celery broker                       |
| `PREDICTION_SERVER`  | `koina.wilhelmlab.org:443` | Koina gRPC endpoint                 |
| `MAX_UPLOAD_BYTES`   | `5000000000`               | Per-file upload limit (5 GB)        |
| `CELERY_CONCURRENCY` | `1`                        | Worker parallelism (1 = FIFO queue) |
| `APP_MODE`           | `local`                    | `local` \| `hosted` (see below)     |

### Custom prediction server

To use an on-prem Koina/Triton instance, set `PREDICTION_SERVER` in `.env`:

```
PREDICTION_SERVER=my-koina-server.local:443
```

Or set it per-job in the **Models & Prediction** form section.

---

## Hardware notes

-   **RAM**: Full-proteome rescoring jobs can require 8–16 GB RAM depending on dataset size.
-   **CPU**: Set `numThreads` in the form (Advanced section) to control parallelism. `CELERY_CONCURRENCY=1` means one job at a time on the worker.
-   **Jobs are long**: Large spectral library generation or rescoring runs can take hours. The status page polls automatically; save the Job ID URL and come back.

---

## Native dependency caveats

### ThermoRawFileParser (`.raw` → `.mzML`)

**Not included** in the Docker image. To use `.raw` Thermo files:

1. Download [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser/releases) and install Mono runtime.
2. Set the `thermoExe` path in the Advanced section of any Rescoring/CE Calibration job (or set a default in `.env`).
3. **Alternative**: pre-convert your `.raw` files to `.mzML` (e.g. with msConvert) and select **mzML** as spectra type — no ThermoRawFileParser needed.

### Percolator

The default FDR method is **mokapot** (pip-installable, included). If you need `percolator` (the binary):

-   Install it from [percolator.ms](http://percolator.ms/) and ensure it's on `$PATH` inside the container.
-   Select `percolator` as the **FDR method** in the Rescoring form.

---

## Running tests

Activate your Python environment, install backend dependencies with Poetry, and then run tests with `poetry run pytest`.

```bash
# Backend unit + integration tests (fast, no Redis needed)
cd backend
poetry run pytest tests/ -v

# End-to-end test (runs a real SpectralLibraryGeneration job)
ENABLE_E2E=1 poetry run pytest tests/test_e2e.py -v -s

# Frontend unit tests
cd frontend && npm run test:unit
```

---

## Finding your Job ID

After submitting, you are redirected to `/jobs/<job_id>`. Bookmark this URL or copy the displayed Job ID. You can also list recent jobs at `/jobs`.

---

## Links

-   [Oktoberfest GitHub](https://github.com/wilhelm-lab/oktoberfest)
-   [Oktoberfest Docs](https://oktoberfest.readthedocs.io)
-   [Koina prediction server](https://koina.wilhelmlab.org)
