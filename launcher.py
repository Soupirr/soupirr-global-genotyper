"""
NDV Genotyper - Entry point for the packaged .exe
"""

import sys
import os
import threading
import time
import webview
from streamlit.web import cli as stcli


def get_base_path():
    """Returns the base directory whether running as exe or as script."""
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    base_path = get_base_path()

    # Set working directory so relative paths (e.g. image/icon.png) resolve correctly
    os.chdir(base_path)

    app_path = os.path.join(base_path, "app.py")

    print("Starting NDV Genotyper...")
    print(f"App path: {app_path}")


def start_streamlit(app_path):
    sys.argv = [
        "streamlit",
        "run",
        app_path,
        "--server.headless=true",
        "--global.developmentMode=false",
        "--server.port=8501",
    ]
    stcli.main()


threading.Thread(target=start_streamlit, args=(app_path,), daemon=True).start()

time.sleep(3)

window = webview.create_window(
    "NDV Genotyper", "http://localhost:8501", width=1400, height=900
)
webview.start()

os._exit(0)
