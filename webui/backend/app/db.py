from __future__ import annotations

from sqlalchemy import create_engine, event, text
from sqlalchemy.orm import DeclarativeBase, sessionmaker

from app.config import settings


def _make_engine():
    url = settings.database_url
    connect_args = {}
    if url.startswith("sqlite"):
        connect_args = {"check_same_thread": False}
    engine = create_engine(url, connect_args=connect_args)
    if url.startswith("sqlite"):
        # Use a high busy timeout to allow concurrent API + worker writes without WAL (NFS support)
        @event.listens_for(engine, "connect")
        def set_sqlite_pragma(conn, _):
            conn.execute("PRAGMA busy_timeout=5000")

    return engine


engine = _make_engine()
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


class Base(DeclarativeBase):
    pass


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
