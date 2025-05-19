"""PyPSA-USA: A Python package for modeling the US energy system using PyPSA."""

# Data handling
from .data import caiso, eia, pudl

# Sector modeling
from .sectors import build_heat, build_transport, electricity

# Visualization
from .visualization import (
    plot_network_maps,
    plot_sankey_carbon,
    plot_sankey_energy,
    plot_statistics,
    summary,
)

# Workflow interface
from .workflow import plot_results, prepare_network, solve_network

__version__ = "0.1.0"  # This should match your project version

__all__ = [
    "build_heat",
    "build_transport",
    "caiso",
    "eia",
    "electricity",
    "plot_network_maps",
    "plot_results",
    "plot_sankey_carbon",
    "plot_sankey_energy",
    "plot_statistics",
    "prepare_network",
    "pudl",
    "solve_network",
    "summary",
]
