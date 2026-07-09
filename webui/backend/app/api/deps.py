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

from typing import Any
from dataclasses import dataclass

from fastapi import Request, Response, HTTPException
import uuid
from dataclasses import dataclass
from app.config import settings

@dataclass
class AppUser:
    id: str
    is_local: bool = True
    is_global: bool = False

LOCAL_USER = AppUser(id="local", is_local=True, is_global=False)

def get_current_user(request: Request, response: Response) -> AppUser:
    """Return the authenticated user.

    Local mode: always returns LocalUser.
    Hosted mode: verify user from request headers (X-User-Id, X-Forwarded-User, Remote-User, or Authorization).
    Fallback to cookie-based UUID sessions.
    """
    if settings.app_mode == "hosted":
        user_id = request.cookies.get("oktoberfest_session")
        if not user_id:
            user_id = str(uuid.uuid4())
            # Set cookie on the response so the browser remembers this UUID
            response.set_cookie(
                key="oktoberfest_session",
                value=user_id,
                httponly=True,
                secure=True,
                max_age=60 * 60 * 24 * 365,  # 1 year
                samesite="lax"
            )

        is_global = False

        # Fallback to checking settings.global_users list
        global_users_list = [u.strip() for u in settings.global_users.split(",") if u.strip()]
        if user_id in global_users_list:
            is_global = True

        return AppUser(id=user_id, is_local=False, is_global=is_global)

    return LOCAL_USER


def check_job_access(job: Any, user: AppUser, read_only: bool = False) -> None:
    """Check if the user has permission to access this job.

    In hosted mode:
    - Global users have access to all jobs.
    - Regular users can access their own jobs.
    - If read_only is True, anyone possessing the job ID can view it.
    """
    if settings.app_mode != "hosted":
        return
    if user.is_global:
        return
    if read_only:
        return
    if job.owner_id != user.id:
        raise HTTPException(status_code=404, detail="Job not found")


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
