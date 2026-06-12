"""
NDV Genotyper - Entry point for the packaged .exe
"""

import sys
import os
import threading
import webbrowser
import time


def get_base_path():
    """Returns the base directory whether running as exe or as script."""
    if getattr(sys, 'frozen', False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


def open_browser():
    """Opens the browser after giving Streamlit time to start."""
    time.sleep(4)
    webbrowser.open('http://localhost:8501')


if __name__ == '__main__':
    base_path = get_base_path()

    # Set working directory so relative paths (e.g. image/icon.png) resolve correctly
    os.chdir(base_path)

    app_path = os.path.join(base_path, 'app.py')

    print("Starting NDV Genotyper...")
    print(f"App path: {app_path}")

    threading.Thread(target=open_browser, daemon=True).start()

    from streamlit.web import cli as stcli
    sys.argv = [
        'streamlit', 'run',
        app_path,
        '--server.headless=true',
        '--global.developmentMode=false',
        '--server.port=8501',
    ]
    sys.exit(stcli.main())
