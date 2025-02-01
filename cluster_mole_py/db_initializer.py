import os
import sqlite3
import pandas as pd
import requests


DB_PATH = "data/project.db"


def create_database():
    """Creates an SQLite database and returns the connection."""
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    print(f"Database created at {DB_PATH}")
    return conn


def load_df_to_db(conn, df, table_name):
    """Loads a dataframe into a database table."""
    df.to_sql(table_name, conn, if_exists="replace", index=False)


def load_panglao_to_db(conn, file_path, table_name):
    df = pd.read_csv(file_path)
    load_df_to_db(conn, df, table_name)
    print("Loaded panglaoDB to table panglaoDB")


def load_cellmarker2_to_db(conn, file_path, table_name):
    df = pd.read_excel(file_path)
    load_df_to_db(conn, df, table_name)
    print("Loaded panglaoDB to table panglaoDB")


def download_pangalodb(save_path):
    url = "https://panglaodb.se/markers.html?cell_type='all_cells'"
    panglaodb_df = pd.read_html(url)[0]
    if panglaodb_df.empty:
        raise Exception("Failed to download PanglaoDB!")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    panglaodb_df.to_csv(save_path, index=False)
    print(f"Downloaded: {save_path}")


def download_cellmarker2(save_path):
    url = "http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_All.xlsx"
    download_files(url, save_path)


def download_files(url, save_path):
    """Downloads a file from a given URL and saves it."""
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    response = requests.get(url, headers)
    if response.status_code == 200:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, "wb") as f:
            f.write(response.content)
        print(f"Downloaded: {save_path}")
    else:
        raise Exception(
            f"Failed to download {url}, Status Code: {response.status_code}"
        )


def setup_database():
    # TODO: Add misigdb for human and mouse
    # TODO: Add scType cell type database
    panglaoDB_file_path = "data/panglaoDB.csv"
    cellmarker2_file_path = "data/cellmarker2.xlsx"

    download_pangalodb(panglaoDB_file_path)
    download_cellmarker2(cellmarker2_file_path)

    conn = create_database()

    load_panglao_to_db(conn, panglaoDB_file_path, "panglaodb")
    load_cellmarker2_to_db(conn, cellmarker2_file_path, "cell_marker")


if __name__ == "__main__":
    setup_database()
