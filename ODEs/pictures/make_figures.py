# generates the figures of the ODE chapter:  python3 make_figures.py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# colors: categorical slots, blue ramp, neutral ink
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
BLUE_LIGHT = "#cde2fb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "savefig.dpi": 200, "savefig.facecolor": "white",
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def axes_cross(ax):
    ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
    ax.axvline(0, color=AXIS, lw=0.8, zorder=1)


# ---------- stability domains of explicit Euler, implicit Euler, Crank-Nicolson

def stability_domains():
    x = np.linspace(-3, 3, 601)
    y = np.linspace(-2.5, 2.5, 501)
    Z = x[None, :] + 1j * y[:, None]
    methods = [
        ("explicit Euler", r"$|1+z| \leq 1$", lambda z: 1 + z, (-1, 0)),
        ("implicit Euler", r"$|1-z| \geq 1$", lambda z: 1 / (1 - z), (-1.6, 1.4)),
        ("Crank-Nicolson", r"$\mathrm{Re}\, z \leq 0$", lambda z: (2 + z) / (2 - z), (-1.5, 1.4)),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(9.5, 3.4), sharey=True)
    for ax, (name, formula, g, pos) in zip(axs, methods):
        with np.errstate(divide="ignore", invalid="ignore"):   # poles on the grid
            G = np.abs(g(Z))
        ax.contourf(x, y, G, levels=[0, 1], colors=[BLUE_LIGHT], zorder=0)
        ax.contour(x, y, G, levels=[1], colors=[BLUE], linewidths=2, zorder=2)
        axes_cross(ax)
        ax.text(*pos, "stable", ha="center", va="center", color=INK, fontsize=10)
        ax.set_title(f"{name}:  {formula}")
        ax.set_aspect("equal")
        ax.set_xlim(-3, 3); ax.set_ylim(-2.5, 2.5)
        ax.set_xlabel(r"$\mathrm{Re}\, z$")
        ax.grid(False)
    axs[0].set_ylabel(r"$\mathrm{Im}\, z$")
    fig.tight_layout()
    fig.savefig("stability_domains.png")
    plt.close(fig)


# ---------- stability regions of explicit Runge-Kutta methods with s = p

def rk_stability():
    x = np.linspace(-3.5, 1.2, 941)
    y = np.linspace(-3.3, 3.3, 1321)
    Z = x[None, :] + 1j * y[:, None]
    colors = [BLUE, ORANGE, AQUA, YELLOW]
    fig, ax = plt.subplots(figsize=(4.8, 5.4))
    axes_cross(ax)
    g = np.ones_like(Z)
    term = np.ones_like(Z)
    for p in range(1, 5):
        term = term * Z / p          # z^p / p!
        g = g + term                 # Taylor polynomial of exp(z)
        cs = ax.contour(x, y, np.abs(g), levels=[1], colors=[colors[p - 1]], linewidths=2, zorder=3)
        # direct label at the top of the curve
        pts = np.vstack([seg for seg in cs.allsegs[0]])
        top = pts[np.argmax(pts[:, 1])]
        ax.annotate(f"$p = {p}$", xy=top, xytext=(top[0] + 0.15, top[1] + 0.12),
                    color=INK, fontsize=10, ha="left", va="bottom")
        ax.plot([], [], color=colors[p - 1], lw=2, label=f"$s = p = {p}$")
    ax.set_aspect("equal")
    ax.set_xlim(-3.5, 1.2); ax.set_ylim(-3.3, 3.3)
    ax.set_xlabel(r"$\mathrm{Re}\, z$"); ax.set_ylabel(r"$\mathrm{Im}\, z$")
    ax.grid(False)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig("rk_stability.png")
    plt.close(fig)


# ---------- explicit Euler for the mass-spring system, as in demos/test_ode.cpp

def massspring():
    tend, steps = 4 * np.pi, 100
    tau = tend / steps
    y = np.zeros((steps + 1, 2))
    y[0] = (1, 0)
    for i in range(steps):
        y[i + 1] = y[i] + tau * np.array([y[i, 1], -y[i, 0]])
    t = np.linspace(0, tend, steps + 1)
    te = np.linspace(0, tend, 1000)

    fig, ax = plt.subplots(figsize=(6, 3.4))
    ax.plot(te, np.cos(te), color=INK2, lw=1.2, ls="--", label="exact solution")
    ax.plot(te, -np.sin(te), color=INK2, lw=1.2, ls="--")
    ax.plot(t, y[:, 0], color=BLUE, lw=2, label="position")
    ax.plot(t, y[:, 1], color=ORANGE, lw=2, label="velocity")
    ax.set_xlabel("time $t$")
    ax.set_xlim(0, tend)
    ax.legend(loc="upper left", ncol=3, fontsize=9, bbox_to_anchor=(0, 1.12))
    fig.tight_layout()
    fig.savefig("massspring_time.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 4.4))
    axes_cross(ax)
    ax.plot(np.cos(te), -np.sin(te), color=INK2, lw=1.2, ls="--", label="exact solution")
    ax.plot(y[:, 0], y[:, 1], color=BLUE, lw=2, label="explicit Euler")
    ax.plot(*y[0], "o", color=BLUE, ms=8, markeredgecolor="white", markeredgewidth=2, zorder=4)
    ax.annotate("start", xy=y[0], xytext=(y[0, 0] + 0.15, y[0, 1] + 0.15), color=INK)
    ax.set_aspect("equal")
    ax.set_xlabel("position"); ax.set_ylabel("velocity")
    ax.legend(loc="lower left", ncol=2, fontsize=9, bbox_to_anchor=(0, 1.0))
    fig.tight_layout()
    fig.savefig("massspring_phase.png")
    plt.close(fig)


# ---------- direction field of y' = f(t,y) with solutions

def rk4_solve(f, t0, y0, T, n=2000):
    t = np.linspace(t0, T, n + 1)
    y = np.empty(n + 1)
    y[0] = y0
    h = t[1] - t[0]
    for i in range(n):
        k1 = f(t[i], y[i]); k2 = f(t[i] + h/2, y[i] + h/2*k1)
        k3 = f(t[i] + h/2, y[i] + h/2*k2); k4 = f(t[i] + h, y[i] + h*k3)
        y[i+1] = y[i] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return t, y


def direction_field():
    f = lambda t, y: np.cos(t) - 0.5*y
    T, ymin, ymax = 10, -2.5, 3
    fig, ax = plt.subplots(figsize=(7, 3.6))
    tt, yy = np.meshgrid(np.linspace(0.25, T - 0.25, 26), np.linspace(ymin + 0.25, ymax - 0.25, 12))
    s = f(tt, yy)
    n = np.sqrt(1 + s**2)
    ax.quiver(tt, yy, 1/n, s/n, color=MUTED, angles="xy", scale=38, width=0.003,
              headwidth=3, headlength=3.5, zorder=1)
    for y0 in (-2, -1, 1, 2.5):
        t, y = rk4_solve(f, 0, y0, T)
        ax.plot(t, y, color=INK2, lw=1, zorder=2)
    t, y = rk4_solve(f, 0, 0, T)
    ax.plot(t, y, color=BLUE, lw=2, zorder=3, label="solution through $(t_0, y_0)$")
    ax.plot(0, 0, "o", color=BLUE, ms=8, markeredgecolor="white", markeredgewidth=2, zorder=4, clip_on=False)
    ax.annotate("$(t_0, y_0)$", xy=(0, 0), xytext=(0.25, -0.55), color=INK,
                bbox=dict(fc="white", ec="none", pad=1), zorder=5)
    ax.plot([], [], color=INK2, lw=1, label="other initial values")
    ax.set_xlim(0, T); ax.set_ylim(ymin, ymax)
    ax.set_xlabel("$t$"); ax.set_ylabel("$y$")
    ax.grid(False)
    ax.legend(loc="lower left", ncol=2, fontsize=9, bbox_to_anchor=(0, 1.0))
    fig.tight_layout()
    fig.savefig("direction_field.png")
    plt.close(fig)


# ---------- explicit Euler for y' = y

def explicit_euler():
    tau, n = 0.5, 4
    t = tau*np.arange(n + 1)
    y = (1 + tau)**np.arange(n + 1)
    te = np.linspace(0, n*tau, 400)
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(te, np.exp(te), color=INK2, lw=1.2, ls="--", label="exact solution $e^t$")
    ax.plot(t, y, color=BLUE, lw=2, label="explicit Euler, $\\tau = 0.5$")
    ax.plot(t, y, "o", color=BLUE, ms=8, markeredgecolor="white", markeredgewidth=2, zorder=4)
    for i in range(n + 1):
        ax.annotate(f"$y_{i}$", xy=(t[i], y[i]), xytext=(t[i] + 0.05, y[i] - 0.55), color=INK)
    # global error at the end
    te_end = n*tau
    ax.annotate("", xy=(te_end + 0.08, np.exp(te_end)), xytext=(te_end + 0.08, y[-1]),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=1))
    ax.text(te_end + 0.14, (np.exp(te_end) + y[-1])/2, "error", color=INK, va="center")
    ax.set_xticks(t, [f"$t_{i}$" for i in range(n + 1)])
    ax.set_xlim(-0.1, n*tau + 0.55); ax.set_ylim(0, 8)
    ax.set_ylabel("$y$")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig("explicit_euler.png")
    plt.close(fig)


# ---------- one step of the improved Euler method for y' = y

def improved_euler():
    tau, y0 = 1.0, 1.0
    ymid = y0 + tau/2*y0          # half step with f(y_n)
    y_ie = y0 + tau*ymid          # full step with f(y_mid)
    y_ee = y0 + tau*y0
    te = np.linspace(0, tau, 200)
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(te, np.exp(te), color=INK2, lw=1.2, ls="--", label="exact solution")
    # the half step to y~ runs along the explicit Euler line (same slope f(y_n))
    ax.plot([0, tau], [y0, y_ee], color=ORANGE, lw=2, label="explicit Euler")
    d = 0.22   # short segment with the slope at the midpoint
    ax.plot([tau/2 - d, tau/2 + d], [ymid - d*ymid, ymid + d*ymid], color=AQUA, lw=2,
            label="slope at the midpoint")
    ax.plot([0, tau], [y0, y_ie], color=BLUE, lw=2, label="improved Euler")
    for (x, yv, c) in [(0, y0, INK2), (tau/2, ymid, AQUA), (tau, y_ie, BLUE), (tau, y_ee, ORANGE)]:
        ax.plot(x, yv, "o", color=c, ms=8, markeredgecolor="white", markeredgewidth=2, zorder=4)
    ax.annotate("$y_n$", xy=(0, y0), xytext=(0.03, y0 - 0.22), color=INK)
    ax.annotate("$\\tilde y$", xy=(tau/2, ymid), xytext=(tau/2 + 0.02, ymid - 0.25), color=INK)
    ax.text(tau/2 + d + 0.03, 1.9, "slope $f(\\tilde y)$", color=INK, va="center")
    ax.annotate("$y_{n+1}$", xy=(tau, y_ie), xytext=(tau + 0.04, y_ie - 0.06), color=INK)
    ax.set_xticks([0, tau/2, tau], ["$t_n$", "$t_n + \\tau/2$", "$t_{n+1}$"])
    ax.set_xlim(-0.05, tau + 0.25); ax.set_ylim(0.7, 3)
    ax.set_ylabel("$y$")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig("improved_euler.png")
    plt.close(fig)


# ---------- one step of the classical RK4 method for y' = -y

def rk4_step():
    f = lambda y: -y
    tau, yn = 1.0, 1.0
    c = [0, tau/2, tau/2, tau]
    ys, ks = [yn], [f(yn)]                       # stage values y^j and slopes k^j
    for j, fac in enumerate([0.5, 0.5, 1.0]):    # y^{j+1} = y_n + fac*tau*k^j
        ys.append(yn + fac*tau*ks[j])
        ks.append(f(ys[-1]))
    slope = (ks[0] + 2*ks[1] + 2*ks[2] + ks[3]) / 6
    ynew = yn + tau*slope
    te = np.linspace(0, tau, 200)

    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    ax.plot(te, np.exp(-te), color=INK2, lw=1.2, ls="--", label="exact solution")
    for j in (1, 2, 3):                          # how the stages are reached from y_n
        ax.plot([0, c[j]], [yn, ys[j]], color=MUTED, lw=1, ls=":", zorder=2)
    ax.plot([], [], color=MUTED, lw=1, ls=":", label="stage $y^j$ from $y_n$ with slope $k^{j-1}$")
    d = 0.11
    for j in range(4):                           # stage slopes as short segments
        ax.plot([c[j] - d, c[j] + d], [ys[j] - d*ks[j], ys[j] + d*ks[j]], color=ORANGE, lw=2, zorder=3)
        ax.plot(c[j], ys[j], "o", color=ORANGE, ms=8, markeredgecolor="white", markeredgewidth=2, zorder=4)
    ax.plot([], [], color=ORANGE, lw=2, marker="o", ms=8, markeredgecolor="white",
            label="stage values $y^j$, slopes $k^j = f(y^j)$")
    ax.plot([0, tau], [yn, ynew], color=BLUE, lw=2, zorder=3, label="RK4 step, mean slope")
    ax.plot(tau, ynew, "o", color=BLUE, ms=8, markeredgecolor="white", markeredgewidth=2, zorder=5)

    ax.annotate("$y_n = y^0$", xy=(0, yn), xytext=(0.03, yn + 0.04), color=INK)
    ax.annotate("$y^1$", xy=(c[1], ys[1]), xytext=(c[1] + 0.03, ys[1] - 0.08), color=INK)
    ax.annotate("$y^2$", xy=(c[2], ys[2]), xytext=(c[2] + 0.03, ys[2] + 0.04), color=INK)
    ax.annotate("$y^3$", xy=(c[3], ys[3]), xytext=(c[3] + 0.04, ys[3] - 0.08), color=INK, va="center")
    ax.annotate("$y_{n+1}$", xy=(tau, ynew), xytext=(tau + 0.05, ynew + 0.01), color=INK, va="center")
    ax.set_xticks([0, tau/2, tau], ["$t_n$", "$t_n + \\tau/2$", "$t_{n+1}$"])
    ax.set_xlim(-0.15, tau + 0.25); ax.set_ylim(0.1, 1.15)
    ax.set_ylabel("$y$")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig("rk4_step.png")
    plt.close(fig)


# ---------- sketch: pendulum with length constraint

def pendulum():
    from matplotlib.patches import Arc, FancyArrowPatch
    l, alpha = 1.0, np.radians(35)
    x0 = np.array([0.0, 0.0])
    x = x0 + l*np.array([np.sin(alpha), -np.cos(alpha)])
    e = (x0 - x) / l                              # unit vector from the mass to the anchor

    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    ax.set_aspect("equal"); ax.axis("off")

    # ceiling with hatching
    ax.plot([-0.45, 0.45], [0, 0], color=INK, lw=1.5)
    for s in np.linspace(-0.42, 0.42, 10):
        ax.plot([s, s + 0.07], [0, 0.07], color=INK2, lw=0.8)

    # admissible positions c(x) = 0
    ax.add_patch(Arc(x0, 2*l, 2*l, theta1=200, theta2=340, color=MUTED, lw=1.2, ls="--"))
    ax.text(-0.98, -0.52, "$c(x) = 0$", color=INK2, ha="right")

    # vertical axis and angle
    ax.plot([0, 0], [0, -1.05*l], color=MUTED, lw=1, ls=":")
    ax.add_patch(Arc(x0, 0.6, 0.6, theta1=270, theta2=270 + np.degrees(alpha), color=INK2, lw=1))
    ax.text(0.39*np.sin(alpha/2), -0.39*np.cos(alpha/2), r"$\alpha$", color=INK, ha="center", va="center")

    # stick and mass
    n = np.array([np.cos(alpha), np.sin(alpha)])  # normal to the stick
    ax.plot([x0[0], x[0]], [x0[1], x[1]], color=INK, lw=2)
    ax.text(*(x0 + 0.3*(x - x0) + 0.08*n), "$l$", color=INK, ha="center", va="center")
    ax.plot(*x0, "o", color=INK, ms=6, zorder=5)
    ax.text(x0[0] - 0.05, x0[1] - 0.12, "$x_0$", color=INK, ha="right")
    ax.plot(*x, "o", color=BLUE, ms=18, markeredgecolor="white", markeredgewidth=2, zorder=5)
    ax.text(x[0] + 0.12, x[1] + 0.02, "$x$, mass $m$", color=INK, va="center")

    # forces
    def arrow(start, vec, color):
        ax.add_patch(FancyArrowPatch(start, start + vec, arrowstyle="-|>", mutation_scale=14,
                                     color=color, lw=2, zorder=4, shrinkA=9, shrinkB=0))
    g = np.array([0, -0.5])
    arrow(x, g, ORANGE)
    ax.text(x[0] + 0.05, x[1] - 0.42, r"$-m g\, e_z$", color=INK)
    arrow(x + 0.07*n, 0.4*e, AQUA)                # parallel to the stick, not on top of it
    ax.text(*(x + 0.07*n + 0.4*e + np.array([0.06, 0.0])), r"$\lambda \nabla_x c$", color=INK, va="center")

    ax.set_xlim(-1.3, 1.35); ax.set_ylim(-1.45, 0.2)
    fig.savefig("pendulum.png", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)


if __name__ == "__main__":
    stability_domains()
    rk_stability()
    massspring()
    direction_field()
    explicit_euler()
    improved_euler()
    rk4_step()
    pendulum()
