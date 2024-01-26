from worker import OKworker
import multiprocessing
from copy import copy
import sys
from utils import GracefulKiller
import time
import os


def spawn(base_dir, sqlite_path):
    """
    Spawns a worker process for the Oktoberfest application.

    :param base_dir: The base directory of the application.
    :param sqlite_path: The path to the SQLite database.
    """
    OKworker(base_dir, sqlite_path)
    sys.exit()


def queue_daemon():
    """
    Function to run the queue daemon.

    This function continuously spawns worker processes to process jobs from a queue.
    It monitors for termination signals and gracefully stops the server when required.
    """
    os.chdir("/home/armin/projects/ok_ui/oktoberfest/oktoberfest/qdaemon")

    sqlite_path = (
        "/home/armin/projects/ok_ui/oktoberfest/oktoberfest/qdaemon/jobs.db"
    )

    killer = GracefulKiller()
    processes = []
    alive = []
    workers = 10
    while not killer.kill_now:
        for i in processes:
            if i.is_alive():
                alive.append(i)
        processes = copy(alive)
        alive = []
        for _i in range(0, workers - len(processes)):
            p = multiprocessing.Process(
                target=spawn, args=("/baseDir", sqlite_path), daemon=True
            )
            p.start()
            processes.append(p)
            time.sleep(0.5)
        time.sleep(30)
    # Process is killed here
    for p in processes:
        if p.is_alive():
            p.terminate()
    print("Server stopped gracefully")


if __name__ == "__main__":
    queue_daemon()
