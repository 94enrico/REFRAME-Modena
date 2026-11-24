"""Reproduce both full manuscript tables, including CD4/CD8 and viral suppression."""

from pathlib import Path
import sys

FOLDER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FOLDER.parent / 'Paper/Scripts'))
from descrittive import write_tables

if __name__ == '__main__':
    write_tables(FOLDER / 'Tables')
