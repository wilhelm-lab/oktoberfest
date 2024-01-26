from worker import OKworker
import multiprocessing
from copy import copy
import sys
from utils import GracefulKiller
import time
import os


def spawn(base_dir, sqlite_path):
    OKworker(base_dir, sqlite_path)
    sys.exit()


def queueDaemon():
    os.chdir("/home/armin/projects/ok_ui/oktoberfest/oktoberfest/qdaemon")

    sqlite_path = "/home/armin/projects/ok_ui/oktoberfest/oktoberfest/qdaemon/jobs.db"

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
        for i in range(0, workers - len(processes)):
            p = multiprocessing.Process(target=spawn, args=("/baseDir", sqlite_path), daemon=True)
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
    queueDaemon()
