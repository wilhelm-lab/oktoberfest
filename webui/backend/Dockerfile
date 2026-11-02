# Build context: the repository root (see webui/docker-compose.yml), so the image runs the oktoberfest
# code of this checkout instead of a PyPI release.

# Stage 1: Build Vue 3 frontend
FROM node:22-alpine AS frontend-build
WORKDIR /app/frontend
COPY webui/frontend/package*.json ./
RUN npm ci --silent
COPY webui/frontend/ ./
RUN npm run build

# Stage 2: Python backend + Oktoberfest
# Pinned to Debian 12 (bookworm): the package names below (libicu72, mono-complete) are release specific
FROM python:3.11-slim-bookworm

# System dependencies
# - mono: ThermoRawFileParser (NOT included; mount it and set thermoExe, or upload mzML)
# - .NET 8 runtime + FlashLFQ: MS1 quantification (quantificationOptions.precursor "flashlfq")
# - git: spectrum_fundamentals / spectrum_io are installed from GitHub
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    git \
    unzip \
    libglib2.0-0 \
    libgomp1 \
    libicu72 \
    mono-complete \
    && rm -rf /var/lib/apt/lists/*

# Install Percolator
RUN curl -fsSL https://github.com/percolator/percolator/releases/download/rel-3-06-01/percolator-v3-06-linux-amd64.deb -o /tmp/percolator.deb && \
    dpkg -i /tmp/percolator.deb && \
    rm /tmp/percolator.deb

# Install the .NET 8 runtime and FlashLFQ
ARG FLASHLFQ_VERSION=2.1.4
RUN curl -fsSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh && \
    bash /tmp/dotnet-install.sh --runtime dotnet --channel 8.0 --install-dir /opt/dotnet && \
    rm /tmp/dotnet-install.sh && \
    curl -fsSL https://github.com/smith-chem-wisc/FlashLFQ/releases/download/${FLASHLFQ_VERSION}/FlashLFQ.zip -o /tmp/FlashLFQ.zip && \
    unzip -q /tmp/FlashLFQ.zip -d /opt/flashlfq && \
    rm /tmp/FlashLFQ.zip
ENV DOTNET_ROOT=/opt/dotnet \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    DOTNET_EXE=/opt/dotnet/dotnet \
    FLASHLFQ_EXE=/opt/flashlfq/CMD.dll

WORKDIR /app

# Upgrade pip first to get reliable JSON index parsing and retry logic
RUN pip install --no-cache-dir --upgrade pip

# Install Python dependencies (split to improve layer caching and retry isolation)
RUN pip install --no-cache-dir --retries 5 --timeout 60 \
    "fastapi>=0.111.0" \
    "uvicorn[standard]>=0.29.0" \
    "sqlalchemy>=2.0.0" \
    "alembic>=1.13.0" \
    "celery[redis]>=5.3.0" \
    "redis>=5.0.0" \
    "pydantic>=2.7.0" \
    "pydantic-settings>=2.3.0" \
    "aiofiles>=23.0.0" \
    "python-multipart>=0.0.9"

# Install oktoberfest from this checkout (with directLFQ for quantification), then move spectrum_fundamentals and
# spectrum_io to SPECTRUM_REF (pyproject pins a feature branch that is already merged into development)
ARG SPECTRUM_REF=development
COPY pyproject.toml README.rst /src/oktoberfest/
COPY oktoberfest/ /src/oktoberfest/oktoberfest/
RUN pip install --no-cache-dir --retries 5 --timeout 120 "/src/oktoberfest[quant]" && \
    pip install --no-cache-dir --no-deps --force-reinstall \
    "spectrum-fundamentals @ git+https://github.com/wilhelm-lab/spectrum_fundamentals.git@${SPECTRUM_REF}" \
    "spectrum-io @ git+https://github.com/wilhelm-lab/spectrum_io.git@${SPECTRUM_REF}"

# Copy backend application code
COPY webui/backend/app/ ./app/
COPY webui/backend/alembic/ ./alembic/
COPY webui/backend/alembic.ini ./

# Copy built frontend into the backend static directory
COPY --from=frontend-build /app/frontend/dist /app/frontend/dist

# Default data dir; writable for any APP_UID so a fresh named volume works
RUN mkdir -p /app/data && chmod 777 /app/data

EXPOSE 8080

# Default: run uvicorn. Worker uses the same image with a different command.
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8080"]
