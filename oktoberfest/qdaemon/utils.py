import sqlite3
import shutil
import signal


def create_database(sqlite_path: str):
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    cursor.execute(
        """
            CREATE TABLE IF NOT EXISTS JOBS
            (
                [ID] INTEGER PRIMARY KEY AUTOINCREMENT,
                [TASK_ID] TEXT UNIQUE,
                [DATE] DATETIME,
                [STATUS] TEXT CHECK(status IN ('UPLOADED', 'PENDING', 'RUNNING', 'DONE', 'FAILED')),
                FOREIGN KEY (ID) REFERENCES CONFIGS(ID)
            )
            """
    )
    cursor.execute(
        """
            CREATE TABLE IF NOT EXISTS CONFIGS
            (
                [ID] INTEGER PRIMARY KEY,
                [CONFIG] BLOB,
                FOREIGN KEY (ID) REFERENCES JOBS(ID)
            )
        """
    )
    conn.commit()
    conn.close()


def print_db_schema(self, sqlite_path: str):
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()

    # Get the table names
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = cursor.fetchall()

    # Print the schema for each table
    for table in tables:
        table_name = table[0]
        print(f"Table: {table_name}")
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = cursor.fetchall()
        for column in columns:
            column_name = column[1]
            column_type = column[2]
            print(f"  {column_name}: {column_type}")

    conn.close()


def zip_folder(folder_path: str, output_path: str):
    shutil.make_archive(base_name=output_path, format="zip", root_dir=folder_path)


class GracefulKiller:
    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, *args, **kwargs):
        self.kill_now = True
