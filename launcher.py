"""
NDV Genotyper - Entry point for the packaged .exe
"""

import sys
import os
import multiprocessing
import time
import urllib.request
import ctypes
import webview


def get_base_path():
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


def streamlit_process(app_path, base_path):
    """Runs Streamlit in its own process so it gets a proper main thread."""
    os.chdir(base_path)
    from streamlit.web import cli as stcli
    sys.argv = [
        "streamlit", "run", app_path,
        "--server.headless=true",
        "--global.developmentMode=false",
        "--server.port=8501",
        "--server.address=127.0.0.1",
    ]
    stcli.main()


def wait_for_streamlit(port=8501, timeout=60):
    start = time.time()
    while time.time() - start < timeout:
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}", timeout=1)
            return True
        except Exception:
            time.sleep(0.5)
    return False


if __name__ == "__main__":
    multiprocessing.freeze_support()  # Required for PyInstaller + multiprocessing

    base_path = get_base_path()
    os.chdir(base_path)
    app_path = os.path.join(base_path, "app.py")

    # Streamlit runs in a separate process (needs its own main thread for signal handlers)
    p = multiprocessing.Process(target=streamlit_process, args=(app_path, base_path))
    p.start()

    wait_for_streamlit()

    def set_icon():
        icon_path = os.path.join(base_path, "image", "icon.png")
        hwnd = ctypes.windll.user32.FindWindowW(None, "NDV Genotyper")
        if hwnd:
            hicon = ctypes.windll.user32.LoadImageW(
                None, icon_path, 1, 0, 0, 0x10  # IMAGE_ICON, LR_LOADFROMFILE
            )
            if hicon:
                ctypes.windll.user32.SendMessageW(hwnd, 0x80, 0, hicon)  # ICON_SMALL
                ctypes.windll.user32.SendMessageW(hwnd, 0x80, 1, hicon)  # ICON_BIG

    # pywebview runs in the main thread of the launcher (required on Windows)
    window = webview.create_window(
        "NDV Genotyper", "http://127.0.0.1:8501", maximized=True
    )
    window.events.loaded += set_icon
    webview.start(gui='edgechromium')

    # Window closed — kill Streamlit and exit cleanly
    p.terminate()
    p.join(timeout=3)
    os._exit(0)
