"""Enable per-panel figure export when VGIC_PANEL_DIR is set (see shared/panel_export.py)."""
import os
import sys
from pathlib import Path

if os.environ.get("VGIC_PANEL_DIR"):
    repo = Path(__file__).resolve().parents[2]
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from shared.panel_export import install_panel_export

    install_panel_export(os.environ["VGIC_PANEL_DIR"], os.environ.get("VGIC_PANEL_LABEL", "figures"))
