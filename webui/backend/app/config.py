from __future__ import annotations

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        # oktoberfest-webui.env in webui/ (docker compose passes it as environment; the second path is for running the
        # backend from webui/backend without docker); environment variables take precedence
        env_file=("../oktoberfest-webui.env", "oktoberfest-webui.env"),
        env_file_encoding="utf-8",
        extra="ignore",
    )

    app_mode: str = "local"  # "local" | "hosted"
    data_dir: str = "./data"
    database_url: str = "sqlite:///./data/app.db"
    redis_url: str = "redis://localhost:6379/0"
    max_upload_bytes: int = 10_000_000_000  # 10 GB
    prediction_server: str = "koina.wilhelmlab.org:443"
    # hosted mode: further Koina servers users may choose (comma-separated host:port); PREDICTION_SERVER is always allowed
    allowed_prediction_servers: str = ""
    cors_origins: str = ""
    # cleanup of old jobs (hourly, app/worker/cleanup.py); 0 disables a rule
    result_retention_days: int = 0
    max_data_bytes: int = 0
    cleanup_min_age_hours: int = 24
    execution_backend: str = "local"
    celery_concurrency: int = 28
    max_threads_per_user: int = 4
    max_total_threads: int = 28
    # hosted mode: unfinished jobs per client (queued, running and drafts of the last 24 h) and upload size per job
    max_open_jobs_per_client: int = 10
    max_job_upload_bytes: int = 100_000_000_000  # 100 GB
    # session ids (cookie oktoberfest_session) that see all jobs in hosted mode; they act as passwords
    global_users: str = ""
    # ThermoRawFileParser for .raw spectra (included in the docker image)
    thermo_exe: str = "ThermoRawFileParser.exe"
    # FlashLFQ command line (CMD.dll) and dotnet for MS1 quantification; empty disables flashlfq
    flashlfq_exe: str = ""
    dotnet_exe: str = "dotnet"


settings = Settings()
