"""Schematic figure: (A) the three copying rules, (B) the four selection scenarios.

Writes the combined figure plus every panel on its own to plots/:
  schematic_full.pdf           (A) + (B), 17.8 cm wide
  schematic_A_strategies.pdf   (A) only
  schematic_B_selection.pdf    (B) only, all four scenarios stacked
  schematic_B<i>_<name>.pdf    one scenario each
PNG copies (300 dpi) are written next to every PDF.

Run:  python3 analysis/illustration/plot_schematic.py
"""
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Polygon, Circle, FancyArrowPatch

OUT = Path(__file__).resolve().parent / "plots"
OUT.mkdir(exist_ok=True)

CM = 1 / 2.54
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Serif",
    "mathtext.it": "Liberation Serif:italic",
    "pdf.fonttype": 42,
})
FS = 8          # text size (pt) in the combined 17.8 cm figure
FS_TAG = 11     # (A) / (B) tags

# ---------------------------------------------------------------- colours
# Phenotype colour scale: low = red/orange, middle = pale, high = dark blue
pheno_cmap = LinearSegmentedColormap.from_list(
    "pheno", ["#d73027", "#f46d43", "#fdae61", "#fee090", "#e0f3f8",
              "#abd9e9", "#74add1", "#4575b4", "#313695"])
# Three demonstrator fish, positions on the phenotype axis (0..1)
POS = {"orange": 0.17, "light": 0.64, "dark": 0.93}
POS_MEAN = np.mean(list(POS.values()))
OPT_POS = 0.55                   # fitness optimum on the phenotype axis
COL = {k: pheno_cmap(v) for k, v in POS.items()}
COL["mean"] = pheno_cmap(POS_MEAN)

ARROW_GREY = "#6e6e6e"
OPT_GREY = "#bdbdbd"
BLUE = "#1f3fe0"
RED = "#e8382f"


def darker(c, f=0.55):
    r, g, b = to_rgb(c)
    return (r * f, g * f, b * f)


# ---------------------------------------------------------------- fish
def _fish_outline():
    """Fish facing right, nose at x = 1, tail fork at x = 0, unit length."""
    # body: tapered ellipse between x = 0.22 (tail joint) and 1 (nose)
    t = np.linspace(0, np.pi, 60)
    bx = 0.61 + 0.39 * np.cos(t)                    # 1 -> 0.22
    taper = 0.55 + 0.45 * (bx - 0.22) / 0.78        # thinner toward the tail
    top = np.column_stack([bx, 0.19 * np.sin(t) * taper ** 0.6])
    bot = top[::-1] * [1, -1]
    tail = np.array([[0.24, 0.02], [0.02, 0.17], [0.09, 0.0],
                     [0.02, -0.17], [0.24, -0.02]])
    return np.vstack([top, tail[:3], tail[3:], bot])


FISH = _fish_outline()
FIN_TOP = np.array([[0.45, 0.15], [0.55, 0.27], [0.68, 0.16]])
FIN_BOT = np.array([[0.48, -0.13], [0.53, -0.21], [0.62, -0.14]])
GILL = np.array([[0.80, 0.11], [0.77, 0.0], [0.80, -0.10]])


def draw_fish(ax, x, y, colour, length=1.0, aspect=1.0, lw=0.7):
    """Draw a fish with its centre at (x, y); `aspect` = y-unit / x-unit
    correction so the fish isn't squashed in non-square data coordinates."""
    def tf(p):
        p = np.asarray(p, float)
        return np.column_stack([x + (p[:, 0] - 0.5) * length,
                                y + p[:, 1] * length * aspect])
    edge = darker(colour)
    for fin in (FIN_TOP, FIN_BOT):
        ax.add_patch(Polygon(tf(fin), closed=True, fc=colour, ec=edge,
                             lw=lw * 0.8, zorder=3, joinstyle="round"))
    ax.add_patch(Polygon(tf(FISH), closed=True, fc=colour, ec=edge, lw=lw,
                         zorder=4, joinstyle="round"))
    ax.plot(*tf(GILL).T, color=edge, lw=lw * 0.7, zorder=5,
            solid_capstyle="round")
    ex, ey = tf([[0.90, 0.045]])[0]
    ax.add_patch(mpl.patches.Ellipse((ex, ey), 0.045 * length,
                                     0.045 * length * aspect, fc="white",
                                     ec=edge, lw=lw * 0.5, zorder=6))


def arrow(ax, x0, x1, y, both=False, lw=2.6, head=9):
    style = f"<|-|>,head_length={head * 0.08},head_width={head * 0.05}" if both \
        else f"-|>,head_length={head * 0.08},head_width={head * 0.05}"
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle=style,
                                 color=ARROW_GREY, lw=lw, zorder=2,
                                 shrinkA=0, shrinkB=0,
                                 mutation_scale=10))


# ---------------------------------------------------------------- panel A
def draw_strategies(ax, fs=FS):
    """Copying rules. Data coords: x 0..10, y 0..13."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 13)
    ax.axis("off")
    # aspect correction so fish keep their shape in this axes
    bb = ax.get_window_extent()
    aspect = (13 / bb.height) / (10 / bb.width)

    rules = [("Copy-the-average", COL["mean"]),
             ("Copy-the-best", COL["light"]),
             ("Copy-the-extreme", COL["orange"])]
    L = 1.35                                         # fish length
    group = [(0.0, 0.6), (-0.85, -0.4), (0.85, -0.4)]  # top, bottom-left, bottom-right
    before = ["light", "orange", "dark"]
    for i, (name, after) in enumerate(rules):
        yc = 12.0 - i * 2.95
        for (dx, dy), key in zip(group, before):
            draw_fish(ax, 1.75 + dx, yc + dy, COL[key], L, aspect)
        for dx, dy in group:
            draw_fish(ax, 8.2 + dx, yc + dy, after, L, aspect)
        arrow(ax, 3.85, 6.1, yc - 0.08)
        ax.text(4.97, yc + 0.95, name, ha="center", va="center", fontsize=fs)

    # fitness optimum over the phenotype colour bar
    x0, x1 = 1.0, 8.9
    ybar0, ybar1 = 1.1, 1.85
    grad = np.linspace(0, 1, 512)[None, :]
    ax.imshow(grad, cmap=pheno_cmap, extent=(x0, x1, ybar0, ybar1),
              aspect="auto", zorder=1, interpolation="bicubic")
    xs = np.linspace(x0, x1, 300)
    mu = x0 + OPT_POS * (x1 - x0)
    ax.plot(xs, 2.15 + 1.75 * np.exp(-0.5 * ((xs - mu) / 0.85) ** 2),
            ls=(0, (4, 3)), color="black", lw=1.2, zorder=2)
    # where the three demonstrators sit on the phenotype axis
    for key in before:
        px = x0 + POS[key] * (x1 - x0)
        ax.plot([px], [ybar1 + 0.12], marker="v", ms=4.2, color=COL[key],
                mec=darker(COL[key]), mew=0.5, zorder=3, clip_on=False)
    ax.text(mu, 4.45, "Fitness Optimum", ha="center", va="center", fontsize=fs)
    ax.text((x0 + x1) / 2, 0.45, "Phenotype Value", ha="center", va="center",
            fontsize=fs)


# ---------------------------------------------------------------- panel B
def gauss(x, mu, sd, h):
    return h * np.exp(-0.5 * ((x - mu) / sd) ** 2)


SCEN = {
    "stabilizing": "Stabilizing\nSelection",
    "directional": "Directional\nSelection",
    "fluctuating": "Fluctuating\nSelection",
    "moving":      "Continuously\nMoving\nOptimum",
}
XL, XR = 2.75, 10.0          # curve range inside a scenario axes
BASE = 0.75                  # baseline height
FW_SD, FW_H = 1.15, 1.55     # fitness function (dashed)
POP_SD = 0.32                # population distribution (solid)


def draw_scenario(ax, kind, fs=FS, label=True):
    """One selection scenario. Data coords: x 0..10, y 0..3."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.2)
    ax.axis("off")
    x = np.linspace(XL, XR, 600)
    lw_fit, lw_pop = 1.3, 1.3
    dash = (0, (4, 2.6))

    def fitness(mu, col, alpha=1.0, sd=FW_SD):
        ax.plot(x, BASE + gauss(x, mu, sd, FW_H), ls=dash, color=col,
                lw=lw_fit, alpha=alpha, zorder=2)

    def opt_line(mu, txt=None, lw=1.2):
        ax.plot([mu, mu], [BASE, BASE + FW_H], color=OPT_GREY, lw=lw, zorder=1)
        if txt:
            ax.text(mu, BASE - 0.36, txt, ha="center", va="center", fontsize=fs)

    def population(mu, sd=POP_SD):
        ax.plot(x, BASE + gauss(x, mu, sd, FW_H), color="black", lw=lw_pop,
                zorder=3)

    if label:
        ax.text(1.3, BASE + 0.75, SCEN[kind], ha="center", va="center",
                fontsize=fs, linespacing=1.15)

    if kind == "stabilizing":
        m = 6.4
        fitness(m, BLUE, sd=1.75)
        population(m, sd=0.6)
        opt_line(m, "Optimum")
        ax.plot([XL, XR], [BASE, BASE], color="black", lw=lw_pop, zorder=3)

    elif kind == "directional":
        m0, m1 = 4.25, 8.55
        fitness(m0, BLUE)
        fitness(m1, RED)
        population(m0, sd=0.24)
        opt_line(m0, "Optimum")
        opt_line(m1, "New Optimum")
        arrow(ax, m0, m1, 2.62)
        ax.text((m0 + m1) / 2, 2.92, r"$\Lambda$", ha="center", va="center",
                fontsize=fs + 1)

    elif kind == "fluctuating":
        m0, m1 = 4.05, 8.75
        fitness(m0, BLUE)
        fitness(m1, RED)
        population((m0 + m1) / 2, sd=0.26)
        opt_line(m0, "Optimum 1")
        opt_line(m1, "Optimum 2")
        arrow(ax, m0, m1, 2.62, both=True)
        ax.text((m0 + m1) / 2, 2.92, r"$\Lambda$, $T$", ha="center",
                va="center", fontsize=fs + 1)

    elif kind == "moving":
        # the optimum slides right at constant speed k; earlier positions
        # fade from blue to the current (red) one, population lags behind
        mus = np.linspace(3.85, 8.75, 4)
        cols = LinearSegmentedColormap.from_list("br", [BLUE, RED])(
            np.linspace(0, 1, len(mus)))
        for i, (m, c) in enumerate(zip(mus, cols)):
            last = i == len(mus) - 1
            fitness(m, c, alpha=1.0 if last or i == 0 else 0.45)
            opt_line(m, lw=1.2 if last or i == 0 else 0.7)
        lag = 1.35
        pop_mu = mus[-1] - lag
        population(pop_mu, sd=0.26)
        ax.text(mus[0], BASE - 0.36, "Optimum(0)", ha="center", va="center",
                fontsize=fs)
        ax.text(mus[-1], BASE - 0.36, "Optimum($t$)", ha="center",
                va="center", fontsize=fs)
        # lag between population mean and current optimum
        ylag = BASE + 0.14
        ax.annotate("", (mus[-1], ylag), (pop_mu, ylag),
                    arrowprops=dict(arrowstyle="<|-|>,head_length=0.3,"
                                    "head_width=0.13", color="#444444",
                                    lw=0.6, shrinkA=0, shrinkB=0), zorder=4)
        ax.text((pop_mu + mus[-1]) / 2 + 0.1, ylag + 0.2, "lag", ha="center",
                va="center", fontsize=fs - 1.5, style="italic",
                color="#444444")
        arrow(ax, mus[0], mus[-1] + 0.35, 2.72)
        ax.text((mus[0] + mus[-1]) / 2, 3.0,
                r"$k$ per generation", ha="center", va="center", fontsize=fs)


# ---------------------------------------------------------------- layouts
def tag(fig, x, y, s, fs=FS_TAG):
    fig.text(x, y, s, ha="left", va="top", fontsize=fs)


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)
    print("wrote", OUT / f"{stem}.pdf")


def full_figure():
    W, H = 17.8 * CM, 11.6 * CM
    fig = plt.figure(figsize=(W, H))
    tag(fig, 0.01, 0.99, "(A)")
    tag(fig, 0.505, 0.99, "(B)")
    axA = fig.add_axes([0.02, 0.0, 0.46, 0.93])
    draw_strategies(axA)
    kinds = list(SCEN)
    h = 0.93 / len(kinds)
    for i, k in enumerate(kinds):
        ax = fig.add_axes([0.52, 0.93 - (i + 1) * h, 0.47, h * 0.97])
        draw_scenario(ax, k)
    save(fig, "schematic_full")


def panel_A():
    fig = plt.figure(figsize=(8.7 * CM, 11.0 * CM))
    ax = fig.add_axes([0.01, 0.0, 0.98, 0.95])
    draw_strategies(ax)
    save(fig, "schematic_A_strategies")


def panel_B():
    fig = plt.figure(figsize=(8.7 * CM, 11.0 * CM))
    h = 1.0 / len(SCEN)
    for i, k in enumerate(SCEN):
        ax = fig.add_axes([0.0, 1.0 - (i + 1) * h, 1.0, h * 0.97])
        draw_scenario(ax, k)
    save(fig, "schematic_B_selection")


def single_scenarios():
    for i, k in enumerate(SCEN, start=1):
        fig = plt.figure(figsize=(8.7 * CM, 2.75 * CM))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        draw_scenario(ax, k)
        save(fig, f"schematic_B{i}_{k}")


if __name__ == "__main__":
    full_figure()
    panel_A()
    panel_B()
    single_scenarios()
