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
    max_upload_bytes: int = 5_000_000_000  # 5 GB
    prediction_server: str = "koina.wilhelmlab.org:443"
    cors_origins: str = ""
    result_retention_days: int = 0
    execution_backend: str = "local"
    celery_concurrency: int = 28
    max_threads_per_user: int = 4
    max_total_threads: int = 28
    # session ids (cookie oktoberfest_session) that see all jobs in hosted mode; they act as passwords
    global_users: str = ""
    # ThermoRawFileParser for .raw spectra (included in the docker image)
    thermo_exe: str = "ThermoRawFileParser.exe"
    # FlashLFQ command line (CMD.dll) and dotnet for MS1 quantification; empty disables flashlfq
    flashlfq_exe: str = ""
    dotnet_exe: str = "dotnet"


settings = Settings()
