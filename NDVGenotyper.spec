# -*- mode: python ; coding: utf-8 -*-
#
# PyInstaller spec file for NDV Genotyper
# Run with: pyinstaller NDVGenotyper.spec
#

from PyInstaller.utils.hooks import collect_all, collect_data_files

# Collect all files from packages that have data files or hidden submodules
st_datas,     st_binaries,     st_hiddenimports     = collect_all('streamlit')
bio_datas,    bio_binaries,    bio_hiddenimports     = collect_all('Bio')
altair_datas, altair_binaries, altair_hiddenimports  = collect_all('altair')
pydeck_datas, pydeck_binaries, pydeck_hiddenimports  = collect_all('pydeck')
pv_datas,     pv_binaries,     pv_hiddenimports      = collect_all('webview')
plotly_datas  = collect_data_files('plotly')

a = Analysis(
    ['launcher.py'],
    pathex=[],
    binaries=(
        st_binaries +
        bio_binaries +
        altair_binaries +
        pydeck_binaries +
        pv_binaries
    ),
    datas=(
        st_datas +
        bio_datas +
        altair_datas +
        pydeck_datas +
        pv_datas +
        plotly_datas +
        [
            # Application source (Streamlit reads app.py as a text file at runtime)
            ('app.py',      '.'),
            ('analyzer.py', '.'),
            # Data files required by the app
            ('data',        'data'),
            ('image',       'image'),
            ('tools',       'tools'),
            # Streamlit theme config
            ('.streamlit',  '.streamlit'),
        ]
    ),
    hiddenimports=(
        st_hiddenimports +
        bio_hiddenimports +
        altair_hiddenimports +
        pydeck_hiddenimports +
        pv_hiddenimports +
        [
            # Plotly
            'plotly',
            'plotly.graph_objects',
            'plotly.express',
            'plotly.subplots',
            # Pandas internals often missed by PyInstaller
            'pandas',
            'pandas._libs.tslibs.np_datetime',
            'pandas._libs.tslibs.nattype',
            'pandas._libs.tslibs.timedeltas',
            'pandas._libs.tslibs.timestamps',
            'pandas._libs.tslibs.offsets',
            # Streamlit extras sometimes missed
            'streamlit.runtime.scriptrunner.magic_funcs',
            'streamlit.components.v1',
            # Tornado (Streamlit's web server)
            'tornado',
            'tornado.web',
            'tornado.websocket',
            'tornado.httpserver',
            # Other
            'webview',
            'webview.platforms.winforms',
            'webview.platforms.edgechromium',
            'pyarrow',
            'click',
            'watchdog',
            'packaging',
            'importlib_metadata',
        ]
    ),
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='NDVGenotyper',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,      # UPX can break some DLLs — keep off for safety
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='image\\icon.png',
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='NDVGenotyper',
)
