import os
import sqlite3
from pathlib import Path
from oktoberfest.runner import run_job
from utils import zip_folder
import time


class OKworker:
    def __init__(self, base_dir, sqlite_path):
        self.pid = os.getpid()
        self.base_dir = Path(base_dir)
        self.sqlite_path = Path(sqlite_path)
        self.conn = sqlite3.connect(self.sqlite_path, isolation_level=None)
        self.conn.execute('PRAGMA journal_mode=wal')
        self.cursor = self.conn.cursor()
        self.check_db()

    def startTransaction(self):
        isDone = False
        while not isDone:
            try:
                # self.conn.execute("PRAGMA locking_mode = RESERVED")
                self.cursor.execute("BEGIN EXCLUSIVE TRANSACTION")
                isDone = True
            except sqlite3.OperationalError:
                print("Lock", self.pid)
                time.sleep(5)
    
    def endTransaction(self):
        self.cursor.execute("COMMIT")
        # self.conn.execute("PRAGMA locking_mode = NORMAL")
        # try:
        #     self.cursor.execute("SELECT COUNT(*) FROM JOBS")
        # except sqlite3.OperationalError:
        #     pass
        self.conn.commit()
        
    # def endTransaction(self):
    #     self.cursor.execute("COMMIT")
    #     self.conn.execute("PRAGMA locking_mode = NORMAL")
    #     self.cursor.execute("")

    def check_db(self):
        # Execute your command here
        self.startTransaction()
        self.cursor.execute("SELECT * FROM JOBS WHERE status='PENDING' ORDER BY ID ASC LIMIT 1;")

        # Fetch the results
        result = self.cursor.fetchone()
        
        if result:
            print(f"Job found: {result}")
            config_path = Path(self.base_dir, result[1], "config.json").resolve()
            if not config_path.exists():
                self.cursor.execute(
                    f"""
                        UPDATE JOBS SET STATUS='FAILED' WHERE ID={result[0]}
                    """
                )
                self.endTransaction()
                self.suicide("ConfigNotFound")
                return -1
            # Update database
            self.cursor.execute(
                f"""
                    UPDATE JOBS SET STATUS='RUNNING' WHERE ID={result[0]}
                """
            )
            self.cursor.execute("COMMIT")
            self.conn.commit()

            # Run Oktoberfest via the function
            if self.run_ok(result, config_path) != 0:
                return -1
        else:
            self.suicide("NoJobFound")
        self.conn.close()
        return 0

    def run_ok(self, result, config_path):
        try:
            run_job(
                config_path=config_path
            )
            outFolder = Path("/output")
            zip_folder(outFolder, outFolder)
            self.startTransaction()
            self.cursor.execute(
                f"""
                    UPDATE JOBS SET STATUS='DONE' WHERE ID={result[0]}
                """
            )
            self.cursor.execute("COMMIT")
            self.conn.commit()
        except Exception as e:
            self.startTransaction()
            self.cursor.execute(
                f"""
                    UPDATE JOBS SET STATUS='FAILED' WHERE ID={result[0]}
                """
            )
            self.cursor.execute("COMMIT")
            self.conn.commit()
            self.suicide(e)
            return -1

    def suicide(self, reason):
        # Cleanup, report, remove PID and suicide
        print(f"Suicide triggered because {reason}")
        self.conn.close()
        
    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    print(os.getcwd())
