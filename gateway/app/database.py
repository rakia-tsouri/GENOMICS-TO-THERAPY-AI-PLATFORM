"""SQLAlchemy engine, session factory, and FastAPI dependency."""
import logging
import time

from sqlalchemy import create_engine, text
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import declarative_base, sessionmaker, Session

from .config import settings

logger = logging.getLogger("gateway.db")

# SQLite needs check_same_thread=False when used with FastAPI's threadpool.
connect_args = (
    {"check_same_thread": False} if settings.database_url.startswith("sqlite") else {}
)

engine = create_engine(settings.database_url, connect_args=connect_args, pool_pre_ping=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()


def get_db():
    """Yield a DB session, closing it after the request."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def wait_for_db(retries: int = 30, delay: float = 2.0) -> None:
    """Block until the database accepts connections.

    Containers can start before Postgres is resolvable/ready; without this the
    gateway would crash on startup. Retries with a fixed backoff instead.
    """
    for attempt in range(1, retries + 1):
        try:
            with engine.connect() as conn:
                conn.execute(text("SELECT 1"))
            if attempt > 1:
                logger.info("Database reachable after %d attempt(s).", attempt)
            return
        except OperationalError as e:
            logger.warning("DB not ready (attempt %d/%d): %s", attempt, retries, str(e).splitlines()[0])
            time.sleep(delay)
    raise RuntimeError(f"Database not reachable after {retries} attempts")
