"""
Universal color palette for colorblind-friendly visualizations.

Based on:
- Wong, B. (2011). Points of view: Color blindness. Nature Methods, 8, 441.
- Okabe & Ito (2002). Color Universal Design (CUD)

These colors are distinguishable for people with:
- Protanopia (red-blind)
- Deuteranopia (green-blind)
- Tritanopia (blue-blind)
"""

# Wong's colorblind-friendly palette
# Reference: https://www.nature.com/articles/nmeth.1618
WONG_PALETTE = {
    'orange': '#E69F00',
    'sky_blue': '#56B4E9',
    'bluish_green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermilion': '#D55E00',
    'reddish_purple': '#CC79A7',
    'black': '#000000',
}

# Ordered list for sequential use
WONG_COLORS = [
    '#0072B2',  # blue
    '#E69F00',  # orange
    '#009E73',  # bluish green
    '#D55E00',  # vermilion
    '#CC79A7',  # reddish purple
    '#56B4E9',  # sky blue
    '#F0E442',  # yellow
]

# Solver-specific colors (using Wong palette)
SOLVER_COLORS = {
    'projection': '#0072B2',  # blue
    'simple': '#D55E00',      # vermilion (instead of red)
    'piso': '#009E73',        # bluish green
}

# TVD limiter markers
LIMITER_MARKERS = {
    'none': 'o',
    'minmod': 's',
    'superbee': '^',
    'vanleer': 'D',
    'mc': 'v',
}

# Line styles for limiters
LIMITER_STYLES = {
    'none': '-',
    'minmod': '--',
    'superbee': '-.',
    'vanleer': ':',
    'mc': (0, (5, 1)),  # loosely dashed
}


def get_color(index: int) -> str:
    """Get color by index (cycles through palette)."""
    return WONG_COLORS[index % len(WONG_COLORS)]


def get_solver_color(solver: str) -> str:
    """Get color for a specific solver."""
    return SOLVER_COLORS.get(solver.lower(), '#000000')


def get_limiter_style(limiter: str) -> tuple:
    """Get marker and line style for a limiter."""
    marker = LIMITER_MARKERS.get(limiter.lower(), 'o')
    linestyle = LIMITER_STYLES.get(limiter.lower(), '-')
    return marker, linestyle


# For matplotlib default cycle
def set_wong_cycle():
    """Configure matplotlib to use Wong's colorblind-friendly palette."""
    import matplotlib.pyplot as plt
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=WONG_COLORS)
