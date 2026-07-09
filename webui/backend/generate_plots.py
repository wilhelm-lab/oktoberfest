import sqlite3
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from datetime import datetime

DB_PATH = "/cmnfs/scratch/oktoberfest_webserver/app.db"
OUTPUT_DIR = "/cmnfs/scratch/oktoberfest_webserver/stats_output"

def main():
    if not os.path.exists(DB_PATH):
        print(f"Error: Database not found at {DB_PATH}")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)

    # Load data into pandas
    query = """
    SELECT 
        id, job_type, status, owner_id, ip_address, 
        created_at, file_count, total_size_bytes, spectra_count
    FROM jobs
    """
    df = pd.read_sql_query(query, conn)
    conn.close()

    if df.empty:
        print("No jobs found in the database.")
        return

    # Convert created_at to datetime
    df['created_at'] = pd.to_datetime(df['created_at'])
    df['date'] = df['created_at'].dt.date
    
    # Fill numeric NaNs with 0
    df['file_count'] = df['file_count'].fillna(0)
    df['total_size_bytes'] = df['total_size_bytes'].fillna(0)
    df['spectra_count'] = df['spectra_count'].fillna(0)

    # Set seaborn style
    sns.set_theme(style="whitegrid")

    # 1. Plot: Cumulative Jobs Over Time
    plt.figure(figsize=(10, 5))
    daily_jobs = df.groupby('date').size().reset_index(name='count')
    daily_jobs['cumulative_count'] = daily_jobs['count'].cumsum()
    sns.lineplot(data=daily_jobs, x='date', y='cumulative_count', marker='o')
    plt.title('Cumulative Jobs Submitted Over Time')
    plt.xlabel('Date')
    plt.ylabel('Total Jobs')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'cumulative_jobs.png'))
    plt.close()

    # 2. Plot: Jobs by Type
    plt.figure(figsize=(10, 5))
    sns.countplot(data=df, x='job_type', order=df['job_type'].value_counts().index, palette="viridis")
    plt.title('Distribution of Jobs by Type')
    plt.xlabel('Job Type')
    plt.ylabel('Count')
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'jobs_by_type.png'))
    plt.close()

    # 3. Plot: Daily Unique Users (by owner_id or IP)
    plt.figure(figsize=(10, 5))
    df['user_identifier'] = df['owner_id'].combine_first(df['ip_address'])
    daily_users = df.groupby('date')['user_identifier'].nunique().reset_index(name='unique_users')
    sns.barplot(data=daily_users, x='date', y='unique_users', color='skyblue')
    plt.title('Daily Unique Users')
    plt.xlabel('Date')
    plt.ylabel('Unique Users')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'daily_unique_users.png'))
    plt.close()

    # 4. Plot: Cumulative Files Processed Over Time
    plt.figure(figsize=(10, 5))
    daily_stats = df.groupby('date')[['file_count', 'spectra_count']].sum().reset_index()
    daily_stats['cum_file_count'] = daily_stats['file_count'].cumsum()
    daily_stats['cum_spectra_count'] = daily_stats['spectra_count'].cumsum()

    sns.lineplot(data=daily_stats, x='date', y='cum_file_count', marker='o', color='tab:blue')
    plt.title('Cumulative Files Processed Over Time')
    plt.xlabel('Date')
    plt.ylabel('Cumulative Files Processed')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'cumulative_files.png'))
    plt.close()

    # 5. Plot: Cumulative Spectra Processed Over Time
    plt.figure(figsize=(10, 5))
    sns.lineplot(data=daily_stats, x='date', y='cum_spectra_count', marker='s', color='tab:orange')
    plt.title('Cumulative Spectra Processed Over Time')
    plt.xlabel('Date')
    plt.ylabel('Cumulative Spectra Processed')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'cumulative_spectra.png'))
    plt.close()

    print(f"Stats generation complete. Plots saved to the '{OUTPUT_DIR}/' directory:")
    for f in os.listdir(OUTPUT_DIR):
        print(f" - {f}")

if __name__ == "__main__":
    main()
