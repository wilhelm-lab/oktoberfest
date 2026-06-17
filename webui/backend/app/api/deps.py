"""
FastAPI dependency seams for multi-user / hosted mode.

In local mode:
  - get_current_user() returns a fixed LocalUser (no auth).
  - rate_limit() is a no-op.
  - check_quota() is a no-op.

In hosted mode (APP_MODE=hosted):
  - Swap get_current_user() to verify an OIDC token / API key.
  - Implement rate_limit() and check_quota() to enforce per-user limits.
  These are placeholders; fill them in without changing any endpoint code.
"""

from __future__ import annotations

from dataclasses import dataclass

from fastapi import Request

from app.config import settings


@dataclass
class AppUser:
    id: str
    is_local: bool = True


LOCAL_USER = AppUser(id="local", is_local=True)


def get_current_user(request: Request) -> AppUser:
    """Return the authenticated user.

    Local mode: always returns LocalUser.
    Hosted mode (TODO): verify Bearer token / session cookie from request.
    """
    if settings.app_mode == "hosted":
        # TODO: implement OIDC/API-key verification here
        # raise HTTPException(401, "Not authenticated") if invalid
        pass
    return LOCAL_USER


async def rate_limit(request: Request) -> None:
    """Per-user/IP rate limiting hook.

    Local mode: no-op.
    Hosted mode (TODO): check rate limit counters.
    """
    if settings.app_mode == "hosted":
        pass  # TODO: implement rate limiting


async def check_quota(user: AppUser) -> None:
    """Per-user job quota check.

    Local mode: no-op.
    Hosted mode (TODO): raise 429 if quota exceeded.
    """
    if settings.app_mode == "hosted":
        pass  # TODO: implement quota check
