"""
NDV Genotyper - Entry point for the packaged .exe
"""

import sys
import os
import ctypes
import multiprocessing
import time
import urllib.request


def get_base_path():
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


def unblock_dlls():
    """Remove Mark of the Web (Zone.Identifier) from bundled DLLs.
    Windows blocks .NET assemblies downloaded from the internet without this."""
    if not getattr(sys, "frozen", False):
        return
    for dirpath, _, filenames in os.walk(sys._MEIPASS):
        for filename in filenames:
            if filename.endswith((".dll", ".pyd")):
                filepath = os.path.join(dirpath, filename)
                ctypes.windll.kernel32.DeleteFileW(f"{filepath}:Zone.Identifier")


def streamlit_process(app_path, base_path):
    """Runs Streamlit in its own process so it gets a proper main thread."""
    os.chdir(base_path)
    from streamlit.web import cli as stcli

    sys.argv = [
        "streamlit",
        "run",
        app_path,
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
    multiprocessing.freeze_support()

    # Unblock DLLs before importing webview (required after downloading from internet)
    unblock_dlls()

    import webview  # Late import — must happen after unblock_dlls()

    base_path = get_base_path()
    os.chdir(base_path)
    app_path = os.path.join(base_path, "app.py")

    p = multiprocessing.Process(target=streamlit_process, args=(app_path, base_path))
    p.start()

    wait_for_streamlit()

    def set_icon():
        icon_path = os.path.join(base_path, "misc", "icon.png")
        hwnd = ctypes.windll.user32.FindWindowW(None, "NDV Genotyper")
        if hwnd:
            hicon = ctypes.windll.user32.LoadImageW(None, icon_path, 1, 0, 0, 0x10)
            if hicon:
                ctypes.windll.user32.SendMessageW(hwnd, 0x80, 0, hicon)
                ctypes.windll.user32.SendMessageW(hwnd, 0x80, 1, hicon)

    window = webview.create_window(
        "NDV Genotyper", "http://127.0.0.1:8501", maximized=True
    )
    window.events.loaded += set_icon
    webview.start(gui="edgechromium")

    p.terminate()
    p.join(timeout=3)
    os._exit(0)
