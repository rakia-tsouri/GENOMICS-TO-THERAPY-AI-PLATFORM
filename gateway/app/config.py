"""Central configuration for the API gateway.

All values are overridable via environment variables (see .env.example at repo root).
"""
from functools import lru_cache
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", extra="ignore")

    # --- App ---
    app_name: str = "Genomics-to-Therapy Gateway"
    api_prefix: str = "/api/v1"
    # Comma-separated list of allowed CORS origins (frontend dev + prod)
    cors_origins: str = "http://localhost:3000,http://127.0.0.1:3000"

    # --- Database ---
    # Defaults to local SQLite so the gateway boots without Postgres for quick dev.
    # docker-compose overrides this with a Postgres URL.
    database_url: str = "sqlite:///./gateway.db"

    # --- Auth / JWT ---
    jwt_secret: str = "change-me-in-production"
    jwt_algorithm: str = "HS256"
    access_token_expire_minutes: int = 60 * 24  # 1 day

    # --- Seed admin (created on first startup if no users exist) ---
    # NOTE: must be a valid, non-reserved email domain (email-validator rejects
    # special-use TLDs like .local/.test), otherwise the seeded admin cannot log in.
    admin_email: str = "admin@platform.io"
    admin_password: str = "admin12345"

    # --- Downstream microservice URLs ---
    genomics_url: str = "http://genomics-validation:8000"
    protein_url: str = "http://protein-analysis-service:8001"
    drug_url: str = "http://drug-discovery-service:8002"
    histopathology_url: str = "http://histopathology-service:8003"
    service_timeout_sec: float = 180.0  # BLAST/ESMFold can be slow

    # --- Uploads (WSI images land on the shared volume the histo service reads) ---
    upload_dir: str = "/data/structures/uploads"
    max_upload_mb: int = 200

    @property
    def cors_origin_list(self) -> list[str]:
        return [o.strip() for o in self.cors_origins.split(",") if o.strip()]


@lru_cache
def get_settings() -> Settings:
    return Settings()


settings = get_settings()
