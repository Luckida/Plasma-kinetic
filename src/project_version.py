"""One release version for the launcher and all output modes."""
from pathlib import Path

VERSION = (Path(__file__).resolve().parent.parent / 'VERSION.txt').read_text(encoding='utf-8').strip()
