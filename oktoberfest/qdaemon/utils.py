import sqlite3
import os
import pandas as pd


def get_db_connection(sqlite_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(sqlite_path)
    conn.row_factory = sqlite3.Row
    return conn


def get_job_list(sqlite_path: str) -> pd.DataFrame:
    conn = get_db_connection(sqlite_path)
    job_list = conn.execute("SELECT * FROM job_list").fetchall()
    conn.close()
    return job_list

def update_job_list(job_list:pd.DataFrame) -> pd.DataFrame:
    for i in job_list.job_id:
        s = job_list.loc[job_list.job_id == i, 'status']
        s_update = update_job(job_id = i, status=s)
        job_list.loc[job_list.job_id == i, 'status'] = s_update
    return job_list

def update_job(job_id:str, status:str) -> str:
    match status:
        case 'waiting':
            

    pass

def start_oktoberfest_job():  # 1a
    pass


def check_running_jobs(data_id:str):  # 1b
    pass


def zip_done_jobs():  # 1b1
    pass


def handle_failed_jobs():  # 1b2
    pass


def start_queue():
    pass
