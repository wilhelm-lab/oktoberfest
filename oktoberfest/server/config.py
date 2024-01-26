import os

WORKDIR = 'workdir'
HASHSIZE = 10
SALT = "SALT"
MAX_CONTENT_LENGTH = 16 * 1000 * 1000
UI_BUILD_DIR = 'oktoberfest/server/dist/'

DB_FILE = os.path.join(WORKDIR, "jobs.db")