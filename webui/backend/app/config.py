from __future__ import annotations

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="ignore")

    app_mode: str = "local"  # "local" | "hosted"
    data_dir: str = "./data"
    database_url: str = "sqlite:///./data/app.db"
    redis_url: str = "redis://localhost:6379/0"
    max_upload_bytes: int = 5_000_000_000  # 5 GB
    prediction_server: str = "koina.wilhelmlab.org:443"
    cors_origins: str = "*"
    result_retention_days: int = 0
    execution_backend: str = "local"
    celery_concurrency: int = 1


settings = Settings()
