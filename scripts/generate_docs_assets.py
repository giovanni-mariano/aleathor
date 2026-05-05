#!/usr/bin/env python3
"""Generate static images used by the documentation."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aleathor as ath


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "assets"


def build_model() -> ath.Model:
    model = ath.Model("docs example")

    sphere = ath.Sphere(0, 0, 0, radius=5.0)
    box = ath.RPP(-10, 10, -10, 10, -10, 10)

    model.add_cell(-sphere, material=1, density=10.5, name="fuel")
    model.add_cell(-box & +sphere, material=2, density=1.0, name="moderator")

    return model


def save_slice_plot(model: ath.Model) -> None:
    _, ax = plt.subplots(figsize=(7, 6))
    ax = model.plot(
        z=0,
        bounds=(-10, 10, -10, 10),
        by_material=True,
        show_colorbar=True,
        contour_by="cell",
        resolution=(300, 300),
        ax=ax,
    )
    ax.set_title("Material slice at z = 0")
    ax.figure.savefig(OUT / "tutorial_slice_materials.png", dpi=160, bbox_inches="tight")
    plt.close(ax.figure)


def save_trace_plot(model: ath.Model) -> None:
    trace = model.trace(start=(-12, 0, 0), end=(12, 0, 0))
    ax = trace.plot(show_legend=True)
    ax.set_title("Ray trace through the model")
    ax.figure.savefig(OUT / "tutorial_trace.png", dpi=160, bbox_inches="tight")
    plt.close(ax.figure)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    model = build_model()
    save_slice_plot(model)
    save_trace_plot(model)


if __name__ == "__main__":
    main()
