# Oktoberfest Web UI

A browser interface for [Oktoberfest](https://github.com/wilhelm-lab/oktoberfest): submit and monitor **Rescoring**,
**Collision Energy Calibration** and **Spectral Library Generation** jobs through guided forms instead of writing JSON
configs by hand.

-   Runs on Windows, macOS and Linux with Docker, on your own machine or as a shared server for a group.
-   Rescoring includes MS1 quantification (search engine or FlashLFQ) and PTM site localization with false
    localization rate (FLR) estimation.
-   Everything runs locally: the Docker image contains Oktoberfest, Percolator, ThermoRawFileParser, FlashLFQ and the
    .NET runtime. Only the peptide property predictions are requested from a [Koina](https://koina.wilhelmlab.org)
    server.

Contents: [Quick start](#quick-start-docker) · [Configuration](#configuration) · [Data](#where-is-my-data) ·
[Shared server](#running-a-shared-server-hosted-mode) · [Included tools](#included-tools) ·
[Without Docker](#running-without-docker-development) · [Tests](#tests) · [Troubleshooting](#troubleshooting)

---

## Quick start (Docker)

### 1. Install the prerequisites

-   **Docker**: Docker Desktop on Windows and macOS, Docker Engine with the Compose plugin (v2.24 or newer) on Linux.
    See [Get Docker](https://docs.docker.com/get-docker/). On Windows, use the WSL2 backend.
-   **Git**: [Download Git](https://git-scm.com/downloads).
-   **Hardware**: at least 16 GB RAM is recommended; full proteome rescoring jobs need 8-16 GB each. The image needs a
    few GB of disk space, plus space for your data.

### 2. Get the code and create your settings file

```bash
git clone https://github.com/wilhelm-lab/oktoberfest.git
cd oktoberfest/webui
cp oktoberfest-webui.env.example oktoberfest-webui.env
mkdir -p data
```

`oktoberfest-webui.env` holds all settings of your installation. Every setting is explained in the file itself.
For a personal installation the defaults work as they are. On Linux, set `APP_UID` and `APP_GID` to the output of
`id -u` and `id -g`, so that the files in `data/` belong to you.

### 3. Start

```bash
docker compose --env-file oktoberfest-webui.env up --build -d
```

The first start builds the image, which takes several minutes. Later starts take seconds.

Always pass `--env-file oktoberfest-webui.env` to `docker compose` (also for `down`, `logs` etc.): the compose file
reads the data folder, user and port from it.

### 4. Open the web UI

Go to **[http://localhost:8080](http://localhost:8080)** (or the `API_PORT` you set).

### 5. Check the installation (optional)

```bash
scripts/smoke_test.sh
```

This checks the API, the worker, the frontend, the quantification and localization options and the tools in the image,
and prints `OK: web UI smoke test passed`. It creates a few empty test jobs.

### Stop, update, logs

```bash
docker compose --env-file oktoberfest-webui.env down                     # stop (jobs and data are kept)
git pull && docker compose --env-file oktoberfest-webui.env up --build -d # update to the latest version
docker compose --env-file oktoberfest-webui.env logs -f worker           # follow the job logs
```

The containers restart automatically after a reboot until you stop them with `down`.

---

## Configuration

All settings are in `oktoberfest-webui.env`; [`oktoberfest-webui.env.example`](oktoberfest-webui.env.example) documents
each of them. Restart after changing a setting: `docker compose --env-file oktoberfest-webui.env up -d`.

| Setting                                    | Default                    | What it does                                                                      |
| ------------------------------------------ | -------------------------- | --------------------------------------------------------------------------------- |
| `HOST_DATA_DIR`                            | `./data`                   | Folder (or Docker volume name) on this machine that stores all jobs               |
| `CONTAINER_DATA_DIR`                       | `/app/data`                | The same folder inside the containers                                             |
| `APP_UID` / `APP_GID`                      | `1000`                     | User and group the containers run as (owner of the job files)                     |
| `API_PORT`                                 | `8080`                     | Port of the web UI                                                                |
| `API_BIND`                                 | `127.0.0.1`                | `127.0.0.1`: only this machine can connect; `0.0.0.0`: all network interfaces     |
| `APP_MODE`                                 | `local`                    | `local` (single user) or `hosted` (shared server, see below)                      |
| `GLOBAL_USERS`                             | empty                      | Hosted mode: admin session ids that see all jobs                                  |
| `CELERY_CONCURRENCY`                       | `1`                        | Jobs that can run at the same time                                                |
| `MAX_TOTAL_THREADS`                        | `4`                        | CPU threads all running jobs may use together                                     |
| `MAX_THREADS_PER_USER`                     | `4`                        | Hosted mode: threads per job                                                      |
| `MAX_OPEN_JOBS_PER_CLIENT`                 | `10`                       | Hosted mode: unfinished jobs per session and client address                       |
| `MAX_JOB_UPLOAD_BYTES`                     | `100000000000`             | Hosted mode: total upload size of one job (100 GB)                                |
| `MAX_UPLOAD_BYTES`                         | `10000000000`              | Maximum size of one uploaded file (10 GB)                                         |
| `RESULT_RETENTION_DAYS`                    | `0` (keep)                 | Delete finished jobs older than this many days                                    |
| `MAX_DATA_BYTES`                           | `0` (no limit)             | When the data folder reaches this size, delete the oldest finished jobs           |
| `CLEANUP_MIN_AGE_HOURS`                    | `24`                       | Jobs younger than this are never deleted by `MAX_DATA_BYTES`                      |
| `PREDICTION_SERVER`                        | `koina.wilhelmlab.org:443` | Koina server pre-filled in the forms; users can change it per job (local mode)    |
| `ALLOWED_PREDICTION_SERVERS`               | empty                      | Hosted mode: further Koina servers users may choose                               |
| `THERMO_EXE`, `FLASHLFQ_EXE`, `DOTNET_EXE` | set by the image           | External tools; only needed [without Docker](#running-without-docker-development) |

To use your own Koina instance, set `PREDICTION_SERVER=my-koina-server:443`. In hosted mode users can only choose
`PREDICTION_SERVER` and the servers in `ALLOWED_PREDICTION_SERVERS`.

**Choosing resources**: each job uses one thread per spectra file (at most 4) and up to 16 GB RAM. For a machine with
32 cores and 128 GB RAM, `CELERY_CONCURRENCY=6` and `MAX_TOTAL_THREADS=24` are a reasonable start. Queued jobs wait
until enough threads are free.

---

## Where is my data?

Everything is stored in `HOST_DATA_DIR` (by default `webui/data/`):

```
data/
  app.db           job database
  .ip_salt         hosted mode only (see below)
  <job_id>/
    inputs/        uploaded files
    output/        Oktoberfest output (see the Oktoberfest documentation, "Generated outputs")
    config.json    the Oktoberfest configuration of the job
    captured.log   log of the run
    results.zip    downloadable archive of the results and config
```

After submitting a job you are redirected to `/jobs/<job_id>`, which updates automatically. Large jobs take hours:
bookmark the page and come back later. `/jobs` lists your recent jobs.

`config.json` can be used to rerun the job with Oktoberfest directly (`python -m oktoberfest --config_path config.json`).
Its paths point into `CONTAINER_DATA_DIR`; set `CONTAINER_DATA_DIR` to the same absolute path as `HOST_DATA_DIR` if they
should also be valid outside the containers.

To back up the web UI, stop it and copy `HOST_DATA_DIR`.

**Deleting old jobs**: by default all jobs are kept. The worker checks every hour:

-   `RESULT_RETENTION_DAYS=30` deletes jobs that finished more than 30 days ago.
-   `MAX_DATA_BYTES=300000000000` deletes, once the data folder reaches 300 GB, finished jobs older than
    `CLEANUP_MIN_AGE_HOURS` (24 h), oldest first, until it is below 300 GB again.

Unsubmitted drafts count as finished from their creation; queued and running jobs are never deleted. Users should
download their results soon after a job finishes.

---

## Running a shared server (hosted mode)

For a server used by several people, set in `oktoberfest-webui.env`:

```env
APP_MODE=hosted
API_BIND=127.0.0.1
HOST_DATA_DIR=/srv/oktoberfest-webui
CELERY_CONCURRENCY=6
MAX_TOTAL_THREADS=24
```

In hosted mode, each browser gets its own session (a cookie) and only lists its own jobs. A job can still be viewed by
anyone who knows its id, but only changed by its owner. Each user can run one job at a time, and the threads per job
are capped by `MAX_THREADS_PER_USER`.

### 1. Create the IP salt

Client IP addresses are only stored as salted hashes (to limit parallel jobs per client). Create the salt once; the
server refuses new jobs without it:

```bash
openssl rand -hex 32 > /srv/oktoberfest-webui/.ip_salt
chmod 600 /srv/oktoberfest-webui/.ip_salt
```

The file must be readable by `APP_UID`. Keep it when you move or back up the data folder.

### 2. Put HTTPS in front of it

The session cookie is only sent over HTTPS (it is marked `Secure`), so without HTTPS browsers cannot keep a session and
jobs cannot be created. Keep `API_BIND=127.0.0.1` and run a reverse proxy with a certificate on the same machine, for
example [Caddy](https://caddyserver.com), which gets certificates automatically:

```
oktoberfest.example.org {
    reverse_proxy 127.0.0.1:8080
}
```

With nginx, also raise the upload limit (`client_max_body_size 10G;`) and pass the client address
(`proxy_set_header X-Forwarded-For $remote_addr;`).

Do not expose port 8080 directly (`API_BIND=0.0.0.0`) on a shared server: the per-client limits trust the
`X-Forwarded-For` header, which only the proxy should set.

Browsers treat `http://localhost` as secure, so for a quick test you can also use an SSH tunnel:
`ssh -L 8080:localhost:8080 your-server`, then open `http://localhost:8080`.

### 3. Admin access (optional)

`GLOBAL_USERS` lists session ids that see and can cancel all jobs. Generate a long random value
(`openssl rand -hex 24`), add it to `GLOBAL_USERS`, restart, and set your browser's `oktoberfest_session` cookie for the
site to that value (developer tools → Application/Storage → Cookies). Anyone who knows the value has admin access.

### 4. Security

What the web UI does to protect a shared server:

-   **Files**: uploads are stored under sanitized names inside the job's folder, and a job can only use its own uploads
    as inputs (paths in a submitted config are ignored). Results and logs are only served from the job's folder.
-   **Network**: only the API port is published, on `127.0.0.1` unless `API_BIND` says otherwise; redis is only
    reachable by the other containers. In hosted mode users cannot choose an arbitrary prediction server, so the server
    cannot be made to connect to other hosts.
-   **Browsers**: hosted mode uses an `HttpOnly`, `Secure`, `SameSite=Lax` session cookie; in local mode only pages
    from `localhost` may call the API from another origin, so other websites cannot read your jobs.
-   **Resources**: hosted mode limits threads per job, running jobs per user, unfinished jobs per client and the upload
    size per file and per job.
-   **Containers**: run as `APP_UID` (redis as its own user) without Linux capabilities and with `no-new-privileges`.

What you need to take care of:

-   Serve it over HTTPS through a reverse proxy and keep `API_BIND=127.0.0.1` (see above). Limit request rates there if
    the server is reachable from the internet.
-   Keep `GLOBAL_USERS` values long and random, and `oktoberfest-webui.env` and `.ip_salt` readable only by you.
-   Anyone with a job id can view that job and download its results (the id is the share link). Do not offer the
    server for confidential data unless access to the site itself is restricted, e.g. by the reverse proxy (VPN, IP
    allowlist or login).
-   Watch the disk space of `HOST_DATA_DIR`, or let old jobs be deleted automatically (`MAX_DATA_BYTES`,
    `RESULT_RETENTION_DAYS`, see [Where is my data?](#where-is-my-data)).
-   Update regularly (`git pull`, rebuild) to get security fixes of the image's packages.

---

## Included tools

The Docker image is built from this checkout, so it always runs the Oktoberfest code of the branch you are on. It
contains:

| Tool                                                                    | Version              | Used for                                   |
| ----------------------------------------------------------------------- | -------------------- | ------------------------------------------ |
| Oktoberfest (with the `quant` extra: directLFQ)                         | this checkout        | all jobs                                   |
| spectrum_fundamentals, spectrum_io                                      | `development` branch | annotation, features, file formats         |
| [Percolator](https://github.com/percolator/percolator)                  | 3.06                 | FDR estimation                             |
| [ThermoRawFileParser](https://github.com/CompOmics/ThermoRawFileParser) | 1.4.4 (with mono)    | `.raw` to mzML conversion                  |
| [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ) + .NET 8        | 2.1.4                | MS1 quantification (`precursor: flashlfq`) |

Versions can be changed with build arguments (`THERMO_RAW_FILE_PARSER_VERSION`, `FLASHLFQ_VERSION`, `SPECTRUM_REF`), e.g.
`docker compose --env-file oktoberfest-webui.env build --build-arg SPECTRUM_REF=main`.

---

## Running without Docker (development)

Prerequisites: Python 3.11, Node 22, Redis, Percolator on `PATH`, and for `.raw` files and FlashLFQ: mono,
ThermoRawFileParser, .NET 8 and FlashLFQ (or upload mzML and use search engine quantification).

```bash
# from the repository root: Oktoberfest from this checkout plus the backend dependencies
python -m venv .venv && source .venv/bin/activate
pip install -e ".[quant]"
pip install fastapi "uvicorn[standard]" sqlalchemy alembic "celery[redis]" redis pydantic pydantic-settings \
    aiofiles python-multipart pytest pytest-asyncio httpx anyio
```

In `webui/oktoberfest-webui.env`, uncomment `DATA_DIR`, `DATABASE_URL` and `REDIS_URL` and set the tool paths
(`THERMO_EXE`, `FLASHLFQ_EXE`, `DOTNET_EXE`). Then start each part in its own terminal:

```bash
redis-server                                             # or: docker run -p 6379:6379 redis:7-alpine
cd webui/backend && uvicorn app.main:app --reload --port 8000
cd webui/backend && celery -A app.worker.celery_app worker -B --loglevel=info  # -B: hourly cleanup
cd webui/frontend && npm ci && npm run dev              # http://localhost:5173 (forwards /api to port 8000)
```

The database tables are created when the API starts.

---

## Tests

```bash
cd webui/backend && python -m pytest tests -q          # backend (a few tests need redis on localhost:6379)
cd webui/backend && ENABLE_E2E=1 python -m pytest tests/test_e2e.py -v  # a real spectral library job (needs Koina)
cd webui/frontend && npx vitest run                    # frontend unit tests
cd webui/frontend && npm run build                     # type check and production build
```

The GitHub workflow `.github/workflows/webui.yml` runs the backend and frontend tests, then builds the Docker image,
starts it in hosted mode with settings from `oktoberfest-webui.env.example` and runs `scripts/smoke_test.sh`.

---

## Troubleshooting

| Problem                                                       | Fix                                                                                                                             |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `env file ... oktoberfest-webui.env not found`                | Create it: `cp oktoberfest-webui.env.example oktoberfest-webui.env`                                                             |
| Settings are ignored, data ends up in a Docker volume         | Pass `--env-file oktoberfest-webui.env` to every `docker compose` command                                                       |
| `Permission denied` / `unable to open database file`          | `HOST_DATA_DIR` must exist and be writable by `APP_UID`/`APP_GID` (`mkdir -p data`, set them to `id -u` / `id -g`)              |
| `port is already allocated`                                   | Another program uses the port: change `API_PORT`                                                                                |
| Hosted mode: jobs cannot be created or disappear after reload | Open the site via HTTPS (or `http://localhost` through an SSH tunnel), see [HTTPS](#2-put-https-in-front-of-it)                 |
| `Job submission is disabled ... missing .ip_salt file`        | [Create the IP salt](#1-create-the-ip-salt) in `HOST_DATA_DIR`                                                                  |
| Other machines cannot connect                                 | `API_BIND=127.0.0.1` only allows this machine: use a reverse proxy (shared server) or `API_BIND=0.0.0.0` (trusted network only) |
| A job fails                                                   | The job page shows its log (`captured.log`); `docker compose --env-file oktoberfest-webui.env logs worker` shows the rest       |

---

## Architecture

```
Browser (Vue 3 + Vuetify)
  │  create job, upload files, submit; poll /api/v1/jobs/{id}
  ▼
api: FastAPI (port 8080, also serves the frontend)
  │  queue (Celery + Redis); the scheduler starts queued jobs when threads are free
  ▼
worker: Celery
  │  writes config.json, runs Oktoberfest, zips output/ into results.zip
  ▼
HOST_DATA_DIR: app.db (SQLite) + one folder per job
```

---

## Links

-   [Oktoberfest on GitHub](https://github.com/wilhelm-lab/oktoberfest)
-   [Oktoberfest documentation](https://oktoberfest.readthedocs.io)
-   [Koina prediction server](https://koina.wilhelmlab.org)
