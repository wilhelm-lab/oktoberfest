"""
Isolate the tests from a local installation: ignore webui/oktoberfest-webui.env and use a temporary data folder.

Runs before the test modules import the app, so the database engine and the data folder are created from these settings.
"""

from __future__ import annotations

import tempfile

import app.config

_defaults = app.config.Settings(_env_file=None)
_data_dir = tempfile.mkdtemp(prefix="oktoberfest-webui-tests-")
_defaults.data_dir = _data_dir
_defaults.database_url = f"sqlite:///{_data_dir}/app.db"

for _name in type(_defaults).model_fields:
    setattr(app.config.settings, _name, getattr(_defaults, _name))
