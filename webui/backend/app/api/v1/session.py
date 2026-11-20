"""
Admin sign-in for hosted mode: a browser becomes an admin (sees and can cancel all jobs) by using one of the
GLOBAL_USERS values as its session id.
"""

from __future__ import annotations

import hmac
import uuid

from fastapi import APIRouter, Depends, HTTPException, Response
from pydantic import BaseModel

from app.api.deps import AppUser, get_current_user, global_users, set_session_cookie
from app.config import settings

router = APIRouter(tags=["session"])


class AdminSignIn(BaseModel):
    key: str


def _session(user: AppUser) -> dict:
    return {"mode": settings.app_mode, "is_admin": user.is_global}


@router.get("")
def get_session(user: AppUser = Depends(get_current_user)):
    """Whether this browser is signed in as admin."""
    return _session(user)


@router.post("/admin")
def admin_sign_in(body: AdminSignIn, response: Response):
    """Make this browser an admin session if the key is one of GLOBAL_USERS."""
    if settings.app_mode != "hosted":
        raise HTTPException(status_code=400, detail="Admin sign-in is only used in hosted mode")
    # compare every key in constant time, so the response time does not reveal how much of a key matched
    matches = [hmac.compare_digest(body.key.encode(), key.encode()) for key in global_users()]
    if not any(matches):
        raise HTTPException(status_code=403, detail="Wrong admin key")
    set_session_cookie(response, body.key)
    return {"mode": settings.app_mode, "is_admin": True}


@router.post("/signout")
def sign_out(response: Response, user: AppUser = Depends(get_current_user)):
    """Leave the admin session: the browser gets a new, normal session (the admin's jobs stay with the admin key)."""
    if settings.app_mode == "hosted" and user.is_global:
        set_session_cookie(response, str(uuid.uuid4()))
    return {"mode": settings.app_mode, "is_admin": False}
