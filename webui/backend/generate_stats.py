import sqlite3
import os

DB_PATH = "/cmnfs/scratch/oktoberfest_webserver/app.db"

def main():
    if not os.path.exists(DB_PATH):
        print(f"Error: Database not found at {DB_PATH}")
        return

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. Total Jobs Submitted
    cursor.execute("SELECT COUNT(*) FROM jobs;")
    total_jobs = cursor.fetchone()[0] or 0

    # 2. Unique Users (using owner_id, or ip_address as fallback if owner_id is not widely used)
    cursor.execute("SELECT COUNT(DISTINCT owner_id) FROM jobs WHERE owner_id IS NOT NULL;")
    unique_owners = cursor.fetchone()[0] or 0
    
    cursor.execute("SELECT COUNT(DISTINCT ip_address) FROM jobs WHERE ip_address IS NOT NULL;")
    unique_ips = cursor.fetchone()[0] or 0

    # 3. Total Files Uploaded
    cursor.execute("SELECT SUM(file_count) FROM jobs;")
    total_files = cursor.fetchone()[0] or 0

    # 4. Total Spectra Processed
    cursor.execute("SELECT SUM(spectra_count) FROM jobs;")
    total_spectra = cursor.fetchone()[0] or 0

    # Print Report
    print("="*40)
    print(" Oktoberfest Usage Statistics ")
    print("="*40)
    print(f"Total Jobs Submitted: {total_jobs}")
    print(f"Total Unique Users (by ID): {unique_owners}")
    print(f"Total Unique IP Addresses: {unique_ips}")
    print(f"Total Files Uploaded: {total_files:,}")
    print(f"Total Spectra Processed: {total_spectra:,}")
    print("="*40)

    conn.close()

if __name__ == "__main__":
    main()
