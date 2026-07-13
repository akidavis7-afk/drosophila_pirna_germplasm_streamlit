from pathlib import Path

def test_streamlit_controls_exist():
    text=Path('streamlit_app.py').read_text(encoding='utf-8')
    assert '.file_uploader(' in text
    assert 'Run piRNA / germ-plasm analysis' in text
    assert 'Download result tables as ZIP' in text
    assert 'SRC_DIR' in text
