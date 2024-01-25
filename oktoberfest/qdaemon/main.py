from utils import periodic_job, GracefulKiller
import time

if __name__ == "__main__":
    sqlite_path = "oktoberfest/results/oktoberfest.sqlite"
    n_workers = 4
    pidfile_path = "oktoberfest/results/oktoberfest.pid"
    killer = GracefulKiller()
    while not killer.kill_now:
        periodic_job(sqlite_path, n_workers, pidfile_path)
        time.sleep(30)
