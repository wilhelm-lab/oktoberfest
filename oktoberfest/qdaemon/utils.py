import sqlite3
import pandas as pd
import shutil
import os
import time
import daemon
import signal
from oktoberfest.runner import run_job

wkdir = "oktoberfest/results"
pidfile_path = os.path.join(wkdir, "oktoberfest.pid")


def periodic_job(sqlite_path: str, n_workers: int, pidfile_path: str):  # 1
    task_list = get_task_list(sqlite_path)
    pending_tasks = check_pending_jobs(task_list)

    n_running_tasks = check_running_jobs(task_list)
    available_workers = n_workers - n_running_tasks

    if pending_tasks is None and available_workers > 0:
        for task in pending_tasks.head(available_workers).iterrows():
            start_oktoberfest_task(id=task.id, pidfile_path=pidfile_path)
    time.sleep(30)


def get_task_list(sqlite_path: str) -> pd.DataFrame:
    conn = get_db_connection(sqlite_path)
    task_list = conn.execute("SELECT * FROM task_list").fetchall()
    conn.close()
    return task_list


def get_db_connection(sqlite_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(sqlite_path)
    conn.row_factory = sqlite3.Row
    return conn


def check_pending_jobs(task_list: pd.DataFrame) -> pd.DataFrame:  # constant
    pending_tasks = task_list.loc[task_list.status == "waiting"]
    if pending_tasks.empty:
        return None
    else:
        return pending_tasks.sort_values(by="id", ascending=True)


def check_running_jobs(task_list: pd.DataFrame) -> int:  # constaant
    return task_list.loc[task_list.status == "running"].dim[0]


def update_task_status_in_db(sqlite_path: str, id: int, new_status: str):
    conn = get_db_connection(sqlite_path)
    conn.execute(
        "UPDATE task_list SET status = ? WHERE id = ?", (new_status, id)
    )
    conn.commit()
    conn.close()


# def start_oktoberfest_task(id: int):
#     p = Process(target=oktoberfest_task, args=(id,))
#     p.start()
#     return p


def start_oktoberfest_task(id: int, pidfile_path: str):
    with daemon.DaemonContext(pidfile=pidfile_path):
        oktoberfest_task(id=id)


def oktoberfest_task(id: int) -> str:  # 1a
    try:
        update_task_status_in_db(id=id, new_status="running")
        run_job(
            config_path=os.path.join(wkdir, "config", str(id) + ".json")
        )  # TODO: change path
        zip_file(id=id)
        update_task_status_in_db(id=id, new_status="done")
        msg = "Task done"
    except:
        update_task_status_in_db(id=id, new_status="failed")
        msg = "Task failed"  # TODO: catch error message

    return msg


def zip_file(id: str):
    shutil.make_archive(
        base_name=str(id) + "result",
        format="zip",
        base_dir=os.path.join(wkdir, "result"),
    )
    pass


class GracefulKiller:
    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self):
        self.kill_now = True
