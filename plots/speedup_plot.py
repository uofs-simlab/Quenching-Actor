import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ---------- helpers ----------

def parse_time(s: str) -> float:
    """Parse '1d-16:38:09' or '20:22:28' into seconds."""
    days = 0
    if "d-" in s:
        d_part, t_part = s.split("d-")
        days = int(d_part)
    else:
        t_part = s
    h, m, sec = map(int, t_part.split(":"))
    return (((days * 24) + h) * 60 + m) * 60 + sec

def amdahl_time(N, T1, s):
    """Amdahl's law for wall time."""
    return T1 * (s + (1.0 - s) / N)

# ---------- raw data ----------

dyn_strings = {
    64:  "1d-16:38:09",
    96:  "1d-03:09:26",
    128: "20:22:28",
    160: "16:21:53",
    256: "10:15:24",
    320: "08:12:18",
}

static64_str = "2d-06:48:38"   # static baseline @ 64 cores

# ---------- convert to numbers ----------

cores = np.array(sorted(dyn_strings.keys()), dtype=float)
dyn_seconds = np.array([parse_time(dyn_strings[int(c)]) for c in cores], dtype=float)
dyn_hours   = dyn_seconds / 3600.0

static64_hours = parse_time(static64_str) / 3600.0

# speedup vs static 64-core baseline
speedup_vs_static = static64_hours / dyn_hours

# ---------- fit Amdahl to dynamic times ----------

p0 = [dyn_hours[0], 0.01]  # initial guess
(params, pcov) = curve_fit(
    amdahl_time,
    cores,
    dyn_hours,
    p0=p0,
    bounds=([0.0, 0.0], [np.inf, 1.0]),
)
T1_hat, s_hat = params
s_err = np.sqrt(np.diag(pcov))[1] if pcov.size else np.nan
Smax = 1.0 / s_hat if s_hat > 0 else np.inf

N_fit = np.linspace(cores.min(), cores.max(), 300)
dyn_fit_hours = amdahl_time(N_fit, T1_hat, s_hat)
speedup_fit_vs_static = static64_hours / dyn_fit_hours

# ideal strong scaling from dynamic 64-core time
ideal_hours = dyn_hours[0] * (cores[0] / cores)

print(f"Fitted serial fraction s = {s_hat:.4e} ± {s_err:.1e}")
print(f"Implied asymptotic speedup 1/s = {Smax:.1f}×")

# ---------- Figure 1: speedup vs static 64 + Amdahl fit ----------

fig1, ax1 = plt.subplots(figsize=(6, 4.5))
ax1.plot(cores, speedup_vs_static, marker="o", label="Measured (dynamic, optimised)")
ax1.plot(N_fit, speedup_fit_vs_static, linestyle="--", label="Amdahl fit")

ax1.set_xlabel("Cores")
ax1.set_ylabel("Speedup vs static 64-core baseline")
ax1.set_title("Speedup with Amdahl Fit")
ax1.grid(True, alpha=0.3)
ax1.legend(loc="upper left")

text = (
    f"s ≈ {s_hat:.2e}\n"
    f"1/s ≈ {Smax:.1f}×"
)
ax1.text(0.04, 0.96, text, transform=ax1.transAxes,
         va="top", ha="left", fontsize=9)

fig1.tight_layout()
fig1.savefig("speedup_vs_static64_amdahl.png", dpi=300)

# ---------- Figure 2: wall-clock time vs cores ----------

fig2, ax2 = plt.subplots(figsize=(6, 4.5))
ax2.plot(cores, dyn_hours, marker="o", label="Measured (dynamic, optimised)")
ax2.plot(cores, ideal_hours, linestyle="--", label="Ideal 1/N from dynamic 64-core")
ax2.plot(N_fit, dyn_fit_hours, linestyle=":", label="Amdahl fit (time)")

ax2.set_xlabel("Cores")
ax2.set_ylabel("Wall-clock time (hours)")
ax2.set_title("Wall-clock Time vs Cores")
ax2.grid(True, alpha=0.3)
ax2.legend(loc="upper right")

fig2.tight_layout()
fig2.savefig("wall_time_vs_cores_amdahl.png", dpi=300)

plt.show()
