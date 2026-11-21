ELEMENTS = ("S1", "S2", "S3", "S4", "linker", "S5", "pore_helix", "filter", "S6_upper", "S6_lower")
SCREEN_SWITCHES = {"wiring": True, "own_place": True, "box_corners": True, "s6_on_own_s5": True,
                   "filter_quarter": False}
class ConfigError(ValueError): pass
class Config: pass
def load_config(*a, **k): raise NotImplementedError
def screen_settings(*a, **k): raise NotImplementedError
