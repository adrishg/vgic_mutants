"""Minimal stand-ins for the names screen.py imports from build_dv.
Only the geometry in screen.py is exercised; nothing here is called."""
STEP1_ERRORS = (OSError,)
STEP2_FILES = ()
class Log:
    def __init__(self): self.lines = []
    def __call__(self, s): self.lines.append(s)
def _git(*a): return ""
def _residues(x): return set()
def _versions(): return {}
def _write_json(*a, **k): raise NotImplementedError
def file_digest(*a, **k): return ""
def lowest_spread_pick(*a, **k): raise NotImplementedError
