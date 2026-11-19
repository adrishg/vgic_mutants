"""Export every subplot of a saved figure as a standalone panel.

Composite figures stay unchanged. When enabled, each ``Figure.savefig`` call
also writes one file per panel so figures can be assembled later without
re-plotting:

* a panel is one axes plus its colorbar axes and any axes stacked on the same
  shared x-axis (for example an RMSF profile with its difference strip);
* figure-level legends and the suptitle are hidden in panels; each
  figure-level legend is written once as ``<stem>__legend<k>``;
* a figure with a single panel is copied as ``<stem>__panel1`` so every
  figure has the same layout in the panel directory.

Enable it for a notebook or script run with::

    VGIC_PANEL_DIR=docs/figures/panels VGIC_PANEL_LABEL=<name> \\
        PYTHONPATH=scripts/panel_export python ...

(``scripts/panel_export/sitecustomize.py`` calls :func:`install_panel_export`.)
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox

_ORIGINAL_SAVEFIG = Figure.savefig


def _slug(text: str) -> str:
    text = re.sub(r"\$[^$]*\$", "", str(text))
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return text[:60]


def _panel_groups(fig: Figure) -> list[list]:
    axes = [ax for ax in fig.axes if ax.get_visible()]
    colorbars = [ax for ax in axes if ax.get_label() == "<colorbar>"]
    mains = [ax for ax in axes if ax not in colorbars]
    groups: list[list] = []
    for ax in mains:
        for group in groups:
            anchor = group[0]
            shared = anchor.get_shared_x_axes().joined(anchor, ax)
            a, b = anchor.get_position(), ax.get_position()
            overlap = min(a.x1, b.x1) - max(a.x0, b.x0)
            # A strip stacked under its parent (shared x, no title of its own),
            # such as an RMSF difference row, belongs to the parent panel;
            # titled axes in a shared-x grid are separate panels.
            if shared and overlap > 0.5 * min(a.width, b.width) and not ax.get_title():
                group.append(ax)
                break
        else:
            groups.append([ax])
    for cbar in colorbars:
        c = cbar.get_position()
        nearest = min(
            groups,
            key=lambda g: min(
                abs(c.x0 - ax.get_position().x1) + abs(c.y0 - ax.get_position().y0)
                for ax in g
            ),
        )
        nearest.append(cbar)
    groups.sort(key=lambda g: (-round(max(ax.get_position().y1 for ax in g), 2),
                               min(ax.get_position().x0 for ax in g)))
    return groups


def _show_shared_axis_labels(group: list):
    """Temporarily restore tick labels and axis labels hidden by shared axes.

    In a ``sharex``/``sharey`` grid only the outer subplots carry tick labels
    and axis titles. A standalone panel needs both, so they are switched on and
    copied from a sibling for the export, then restored.
    """
    undo = []
    for ax in group:
        if ax.get_label() == "<colorbar>":
            continue
        for axis, side, getter, setter, shared in (
            (ax.yaxis, "labelleft", ax.get_ylabel, ax.set_ylabel, ax.get_shared_y_axes()),
            (ax.xaxis, "labelbottom", ax.get_xlabel, ax.set_xlabel, ax.get_shared_x_axes()),
        ):
            siblings = [s for s in shared.get_siblings(ax) if s is not ax]
            if not siblings:
                continue
            key = "label1On"
            was_on = axis._major_tick_kw.get(key, True)
            if not was_on:
                ax.tick_params(axis="y" if axis is ax.yaxis else "x", **{side: True})
                undo.append(lambda ax=ax, axis=axis, side=side: ax.tick_params(
                    axis="y" if axis is ax.yaxis else "x", **{side: False}))
            if not getter():
                label = next((s.get_ylabel() if axis is ax.yaxis else s.get_xlabel()
                              for s in siblings
                              if (s.get_ylabel() if axis is ax.yaxis else s.get_xlabel())), "")
                if label:
                    setter(label)
                    undo.append(lambda setter=setter: setter(""))
    def restore():
        for action in reversed(undo):
            action()
    return restore


def export_panels(fig: Figure, out_dir: str | Path, stem: str, dpi: int = 600) -> list[Path]:
    """Write each panel (and figure legend) of ``fig`` to ``out_dir``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    renderer = fig.canvas.get_renderer()
    groups = _panel_groups(fig)
    figure_artists = list(fig.legends) + list(fig.texts)
    if fig._suptitle is not None:
        figure_artists.append(fig._suptitle)
    written = []
    all_axes = [ax for g in groups for ax in g]
    for index, group in enumerate(groups, 1):
        title = next((ax.get_title() for ax in group if ax.get_title()), "")
        name = f"{stem}__panel{index}" + (f"_{_slug(title)}" if title else "")
        hidden = [ax for ax in all_axes if ax not in group] + figure_artists
        states = [artist.get_visible() for artist in hidden]
        for artist in hidden:
            artist.set_visible(False)
        restore_axes = _show_shared_axis_labels(group)
        try:
            box = Bbox.union([ax.get_tightbbox(renderer) for ax in group])
            box = box.transformed(fig.dpi_scale_trans.inverted()).padded(0.06)
            for ext in ("png", "pdf"):
                path = out_dir / f"{name}.{ext}"
                _ORIGINAL_SAVEFIG(fig, path, dpi=dpi, bbox_inches=box, facecolor="white")
                written.append(path)
        finally:
            for artist, state in zip(hidden, states):
                artist.set_visible(state)
            restore_axes()
    for index, legend in enumerate(fig.legends, 1):
        box = legend.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted()).padded(0.06)
        hidden = all_axes + [a for a in figure_artists if a is not legend]
        states = [artist.get_visible() for artist in hidden]
        for artist in hidden:
            artist.set_visible(False)
        try:
            path = out_dir / f"{stem}__legend{index}.png"
            _ORIGINAL_SAVEFIG(fig, path, dpi=dpi, bbox_inches=box, facecolor="white")
            written.append(path)
        finally:
            for artist, state in zip(hidden, states):
                artist.set_visible(state)
    return written


def install_panel_export(out_root: str | Path, label: str) -> Path:
    """Patch ``Figure.savefig`` so every saved figure also writes its panels."""
    target = Path(out_root) / label
    if getattr(Figure.savefig, "_vgic_panel_hook", False):
        return target
    seen: set[tuple[int, str]] = set()

    def savefig(self, fname, *args, **kwargs):
        result = _ORIGINAL_SAVEFIG(self, fname, *args, **kwargs)
        try:
            stem = Path(str(fname)).stem
        except TypeError:  # file-like objects (e.g. inline display buffers)
            return result
        key = (id(self), stem)
        if key not in seen and str(fname).lower().endswith((".png", ".pdf", ".svg")):
            seen.add(key)
            try:
                export_panels(self, target, stem)
            except Exception as error:  # never break the analysis because of panel export
                print(f"[panel_export] {stem}: {type(error).__name__}: {error}")
        return result

    savefig._vgic_panel_hook = True
    Figure.savefig = savefig
    return target
