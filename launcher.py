"""
NDV Genotyper - Entry point for the packaged .exe
"""

import sys
import os
import threading
import time
import urllib.request
import webview
from streamlit.web import cli as stcli


def get_base_path():
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


def start_streamlit(app_path):
    sys.argv = [
        "streamlit", "run", app_path,
        "--server.headless=true",
        "--global.developmentMode=false",
        "--server.port=8501",
    ]
    stcli.main()


def wait_for_streamlit(port=8501, timeout=60):
    """Polls until Streamlit responds or timeout is reached."""
    start = time.time()
    while time.time() - start < timeout:
        try:
            urllib.request.urlopen(f"http://localhost:{port}", timeout=1)
            return True
        except Exception:
            time.sleep(0.5)
    return False


if __name__ == "__main__":
    base_path = get_base_path()
    os.chdir(base_path)
    app_path = os.path.join(base_path, "app.py")

    threading.Thread(target=start_streamlit, args=(app_path,), daemon=True).start()

    wait_for_streamlit()

    window = webview.create_window(
        "NDV Genotyper", "http://localhost:8501", width=1400, height=900
    )
    webview.start()

    os._exit(0)
