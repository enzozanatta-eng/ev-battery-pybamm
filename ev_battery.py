"""
EV battery cell: electrochemical modelling with PyBaMM
=======================================================
Cell: LG M50 21700 (5 Ah, NMC811 / graphite-SiOx), parameter set "Chen2020"

  1. Model hierarchy: Single Particle Model (SPM), SPM with electrolyte (SPMe)
     and the full Doyle-Fuller-Newman model (DFN), compared at 0.5C to 3C
  2. Rate capability and heating: delivered capacity, energy and temperature
     rise for three cooling levels, heat source breakdown and an energy-balance check
  3. Electric vehicle application: a simple vehicle model turns speed into pack
     power and cell power; range at constant highway speeds, the effect of
     aerodynamic drag on range, and a mixed city/highway cycle

Units: SI unless noted. Requires: pip install pybamm matplotlib
"""
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pybamm

pybamm.set_logging_level("ERROR")

# ============================ SETTINGS ============================
PARAMS = "Chen2020"
C_RATES = [0.5, 1, 1.5, 2]
COOLING = {"Natural convection": 10.0, "Forced air": 40.0, "Liquid cooling": 150.0}   # h [W/m2.K]

# Vehicle (mid-size electric car) and pack
MASS = 1900.0          # kg, with driver
CD, AREA = 0.28, 2.3   # drag coefficient, frontal area [m2]
CRR = 0.009            # rolling resistance coefficient
RHO_AIR, G = 1.2, 9.81
ETA_DRIVE = 0.90       # battery-to-wheel efficiency
ETA_REGEN = 0.60       # wheel-to-battery efficiency when braking
P_AUX = 500.0          # auxiliaries [W]
N_SERIES, N_PARALLEL = 96, 46      # 4416 cells, about 400 V and 80 kWh
SPEEDS = [80, 100, 120, 140]       # km/h
OUT = "results"
# ==================================================================

os.makedirs(OUT, exist_ok=True)
N_CELLS = N_SERIES * N_PARALLEL


def params(h=None):
    p = pybamm.ParameterValues(PARAMS)
    if h is not None:
        p["Total heat transfer coefficient [W.m-2.K-1]"] = h
    return p


def discharge(model, c_rate, p):
    exp = pybamm.Experiment([f"Discharge at {c_rate}C until 2.5 V"])
    t0 = time.perf_counter()
    sol = pybamm.Simulation(model, parameter_values=p, experiment=exp).solve()
    return sol, time.perf_counter() - t0


def wheel_power(v, a, cd=CD):
    """Battery power [W] for speed v [m/s] and acceleration a [m/s2]."""
    f = MASS * a + 0.5 * RHO_AIR * cd * AREA * v ** 2 + CRR * MASS * G * (v > 0.1)
    p_wheel = f * v
    return np.where(p_wheel >= 0, p_wheel / ETA_DRIVE, p_wheel * ETA_REGEN) + P_AUX


def mixed_cycle():
    """A 30-minute synthetic cycle: city stop-and-go, suburban and highway (speed in km/h, 1 s steps)."""
    pts = [(0, 0), (20, 50), (60, 50), (75, 0), (95, 0), (115, 40), (170, 40), (185, 0), (205, 0),
           (230, 60), (330, 60), (350, 30), (420, 30), (440, 0), (460, 0), (500, 80), (650, 80),
           (700, 110), (1000, 110), (1100, 120), (1450, 120), (1520, 70), (1650, 70), (1700, 30),
           (1760, 30), (1800, 0)]
    t = np.arange(0, 1801.0)
    v = np.interp(t, [p[0] for p in pts], [p[1] for p in pts])
    return t, v


def run_power_profile(p_cell, dt, model, p, label):
    """Repeat a cell power profile until the cell reaches 2.5 V."""
    t = np.arange(len(p_cell)) * dt
    data = np.c_[t, p_cell]
    exp = pybamm.Experiment([pybamm.step.power(data, termination="< 2.5 V", period=dt)])
    return pybamm.Simulation(model, parameter_values=p, experiment=exp).solve()


# ==================================================================
def main():
    report = []

    def log(text=""):
        print(text)
        report.append(text)

    p0 = params()
    cap_nom = p0["Nominal cell capacity [A.h]"]
    log(f"CELL: LG M50 21700, {cap_nom:.0f} Ah, parameter set {PARAMS}")
    log(f"PACK: {N_SERIES}s{N_PARALLEL}p = {N_CELLS} cells\n")

    # ---------------- 1. Model hierarchy ----------------
    models = {"SPM": pybamm.lithium_ion.SPM, "SPMe": pybamm.lithium_ion.SPMe, "DFN": pybamm.lithium_ion.DFN}
    hier = {}
    for c in C_RATES:
        for name, cls in models.items():
            sol, dt_run = discharge(cls(), c, params())
            hier[(name, c)] = dict(q=sol["Discharge capacity [A.h]"].entries, v=sol["Voltage [V]"].entries,
                                   cap=sol["Discharge capacity [A.h]"].entries[-1], time=dt_run)
    log("MODEL HIERARCHY (isothermal 25 C): voltage error vs the DFN model")
    log(f"{'C-rate':>6} {'model':>5} {'capacity [Ah]':>14} {'RMS error [mV]':>15} {'run time [s]':>13}")
    for c in C_RATES:
        ref = hier[("DFN", c)]
        for name in models:
            r = hier[(name, c)]
            qmax = 0.90 * min(r["q"][-1], ref["q"][-1])     # skip the steep knee before the cut-off
            qq = np.linspace(0, qmax, 300)
            err = np.sqrt(np.mean((np.interp(qq, r["q"], r["v"]) - np.interp(qq, ref["q"], ref["v"])) ** 2)) * 1e3
            r["rms"] = err
            log(f"{c:>6} {name:>5} {r['cap']:>14.3f} {err:>15.1f} {r['time']:>13.2f}")

    # ---------------- 2. Rate capability and heating ----------------
    therm = {}
    for cname, h in COOLING.items():
        for c in C_RATES:
            sol, _ = discharge(pybamm.lithium_ion.DFN(options={"thermal": "lumped"}), c, params(h))
            t = sol["Time [s]"].entries
            T = sol["X-averaged cell temperature [K]"].entries - 273.15
            I = sol["Current [A]"].entries
            V = sol["Voltage [V]"].entries
            therm[(cname, c)] = dict(sol=sol, t=t, T=T, cap=sol["Discharge capacity [A.h]"].entries[-1],
                                     energy=np.trapezoid(I * V, t) / 3600, Tmax=T.max())
    log("\nRATE CAPABILITY AND HEATING (DFN with lumped thermal model, ambient 25 C)")
    log(f"{'cooling':<20} {'C-rate':>6} {'capacity [Ah]':>14} {'energy [Wh]':>12} {'T max [C]':>10}")
    for cname in COOLING:
        for c in C_RATES:
            r = therm[(cname, c)]
            log(f"{cname:<20} {c:>6} {r['cap']:>14.3f} {r['energy']:>12.2f} {r['Tmax']:>10.1f}")

    # energy balance check at 2C, natural convection
    s = therm[("Natural convection", 2)]["sol"]
    t = s["Time [s]"].entries
    I, V = s["Current [A]"].entries, s["Voltage [V]"].entries
    q_irr = s["Ohmic heating [W]"].entries + s["Irreversible electrochemical heating [W]"].entries
    e_heat = np.trapezoid(q_irr, t)
    e_surf = np.trapezoid(I * (s["Surface open-circuit voltage [V]"].entries - V), t)
    e_bulk = np.trapezoid(I * (s["Bulk open-circuit voltage [V]"].entries - V), t)
    log("\nENERGY BALANCE CHECK (2C, natural convection)")
    log(f"  electrical loss I(OCV_surface - V) = {e_surf:.0f} J vs ohmic + irreversible heat = {e_heat:.0f} J "
        f"({(e_heat / e_surf - 1) * 100:+.1f}%)")
    log(f"  with the bulk OCV the loss is {e_bulk:.0f} J: the extra {e_bulk - e_surf:.0f} J is energy stored in the "
        "lithium concentration gradients inside the particles, not yet released as heat")
    q_cap = np.trapezoid(I, t) / 3600
    log(f"  Coulomb counting: integral of current = {q_cap:.4f} Ah vs reported capacity "
        f"{s['Discharge capacity [A.h]'].entries[-1]:.4f} Ah")

    # ---------------- 3. EV application ----------------
    log("\nELECTRIC VEHICLE: range at constant speed (SPMe with lumped thermal, natural convection)")
    log(f"  vehicle {MASS:.0f} kg, Cd {CD}, A {AREA} m2, Crr {CRR}, pack {N_SERIES}s{N_PARALLEL}p")
    ev = {}
    for cd in (CD, CD * 0.9):
        for kmh in SPEEDS:
            v = kmh / 3.6
            p_pack = float(wheel_power(np.array([v]), np.array([0.0]), cd)[0])
            p_cell = p_pack / N_CELLS
            sol = run_power_profile(np.full(3600, p_cell), 10.0, pybamm.lithium_ion.SPMe(options={"thermal": "lumped"}),
                                    params(COOLING["Natural convection"]), f"{kmh}")
            t_end = sol["Time [s]"].entries[-1]
            ev[(cd, kmh)] = dict(range_km=v * t_end / 1000, p_pack=p_pack, p_cell=p_cell, hours=t_end / 3600,
                                 Tmax=sol["X-averaged cell temperature [K]"].entries.max() - 273.15)
    log(f"{'speed [km/h]':>12} {'pack power [kW]':>16} {'cell power [W]':>15} {'range [km]':>11} "
        f"{'range Cd -10%':>14} {'gain':>6}")
    for kmh in SPEEDS:
        a, b = ev[(CD, kmh)], ev[(CD * 0.9, kmh)]
        log(f"{kmh:>12} {a['p_pack'] / 1e3:>16.1f} {a['p_cell']:>15.2f} {a['range_km']:>11.0f} "
            f"{b['range_km']:>14.0f} {(b['range_km'] / a['range_km'] - 1) * 100:>5.1f}%")

    tc, vc = mixed_cycle()
    vms = vc / 3.6
    acc = np.gradient(vms, tc)
    p_pack_cycle = wheel_power(vms, acc)
    n_rep = 12
    p_cell_cycle = np.tile(p_pack_cycle / N_CELLS, n_rep)
    mix = run_power_profile(p_cell_cycle, 1.0, pybamm.lithium_ion.SPMe(options={"thermal": "lumped"}),
                            params(COOLING["Natural convection"]), "mixed")
    t_mix = mix["Time [s]"].entries
    dist = np.cumsum(np.tile(vms, n_rep))[: len(np.arange(int(t_mix[-1])))] / 1000
    range_mix = float(np.trapezoid(np.interp(np.arange(t_mix[-1]), np.arange(len(vms) * n_rep), np.tile(vms, n_rep)))) / 1000
    cycle_km = np.trapezoid(vms, tc) / 1000
    log(f"\nMIXED CYCLE: {cycle_km:.1f} km per 30-minute cycle, average {cycle_km / 0.5:.0f} km/h, "
        f"peak pack power {p_pack_cycle.max() / 1e3:.0f} kW, peak regen {p_pack_cycle.min() / 1e3:.0f} kW")
    log(f"  range {range_mix:.0f} km ({t_mix[-1] / 3600:.1f} h), max cell temperature "
        f"{mix['X-averaged cell temperature [K]'].entries.max() - 273.15:.1f} C")

    # ================= FIGURES =================
    # 1. model hierarchy
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    styles = {"SPM": ":", "SPMe": "--", "DFN": "-"}
    cols = plt.get_cmap("plasma")(np.linspace(0.05, 0.85, len(C_RATES)))
    for c, col in zip(C_RATES, cols):
        for name in models:
            r = hier[(name, c)]
            ax.plot(r["q"], r["v"], styles[name], color=col, lw=1.6,
                    label=f"{name}, {c}C" if c == C_RATES[-1] or name == "DFN" else None)
    ax.set_xlabel("Discharge capacity [Ah]")
    ax.set_ylabel("Voltage [V]")
    ax.set_title("Discharge curves (solid DFN, dashed SPMe, dotted SPM)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    ax = axes[1]
    for name, mk in zip(models, "os^"):
        if name == "DFN":
            continue
        ax.plot(C_RATES, [hier[(name, c)]["rms"] for c in C_RATES], mk + "-", label=name)
    ax.set_xlabel("C-rate")
    ax.set_ylabel("RMS voltage error vs DFN [mV]")
    ax.set_title("Simplified models lose accuracy at high current")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUT}/model_hierarchy.png", dpi=150)
    plt.close(fig)

    # 2. temperature and capacity
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    for cname, ls in zip(COOLING, ["-", "--", ":"]):
        for c, col in zip(C_RATES, cols):
            r = therm[(cname, c)]
            ax.plot(r["t"] / 60, r["T"], ls, color=col, label=f"{c}C" if cname == "Natural convection" else None)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Cell temperature [C]")
    ax.set_title("Temperature during discharge (solid: natural, dashed: air, dotted: liquid)")
    ax.title.set_fontsize(9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax = axes[1]
    for cname, mk in zip(COOLING, "os^"):
        ax.plot(C_RATES, [therm[(cname, c)]["Tmax"] for c in C_RATES], mk + "-", label=cname)
    ax.axhline(60, color="tab:red", ls="--", lw=1)
    ax.text(C_RATES[0], 61, "60 C: typical upper limit for discharge", color="tab:red", fontsize=8)
    ax.set_xlabel("C-rate")
    ax.set_ylabel("Peak cell temperature [C]")
    ax.set_title("Cooling sets how hard the cell can be pushed")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/thermal.png", dpi=150)
    plt.close(fig)

    # 3. heat sources at 2C
    fig, ax = plt.subplots(figsize=(8, 4.5))
    tt = s["Time [s]"].entries / 60
    srcs = [("Ohmic heating [W]", "Ohmic (electronic + ionic resistance)"),
            ("Irreversible electrochemical heating [W]", "Irreversible (reaction overpotential)")]
    for key, lab in srcs:
        ax.plot(tt, s[key].entries, label=lab)
    ax.plot(tt, s["Total heating [W]"].entries, "k", lw=2, label="Total")
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Heat generation [W]")
    ax.set_title("Where the heat comes from (2C discharge)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/heat_sources.png", dpi=150)
    plt.close(fig)

    # 4. EV range
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    ax.plot(SPEEDS, [ev[(CD, k)]["range_km"] for k in SPEEDS], "o-", color="k", label=f"Cd = {CD}")
    ax.plot(SPEEDS, [ev[(CD * 0.9, k)]["range_km"] for k in SPEEDS], "s--", color="tab:green",
            label=f"Cd = {CD * 0.9:.3f} (10% less drag)")
    ax.set_xlabel("Constant speed [km/h]")
    ax.set_ylabel("Range [km]")
    ax.set_title("Range falls fast with speed: aerodynamic drag grows with V squared")
    ax.title.set_fontsize(10)
    ax.grid(alpha=0.3)
    ax.legend()
    ax = axes[1]
    ax.plot(tc / 60, vc, color="tab:blue")
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Speed [km/h]", color="tab:blue")
    ax2 = ax.twinx()
    ax2.plot(tc / 60, p_pack_cycle / 1e3, color="tab:orange", lw=0.9)
    ax2.set_ylabel("Pack power [kW] (negative = regenerative braking)", color="tab:orange")
    ax.set_title(f"Mixed cycle: {cycle_km:.1f} km in 30 min")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/ev_range.png", dpi=150)
    plt.close(fig)

    # 5. mixed cycle battery response
    fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
    th = t_mix / 3600
    axes[0].plot(th, mix["Voltage [V]"].entries, lw=0.7)
    axes[0].set_ylabel("Cell voltage [V]")
    axes[1].plot(th, mix["Current [A]"].entries, lw=0.5, color="tab:orange")
    axes[1].set_ylabel("Cell current [A]")
    axes[2].plot(th, mix["X-averaged cell temperature [K]"].entries - 273.15, color="tab:red")
    axes[2].set_ylabel("Temperature [C]")
    axes[2].set_xlabel("Time [h]")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle(f"Cell response over repeated mixed cycles until empty: range {range_mix:.0f} km", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{OUT}/mixed_cycle_response.png", dpi=150)
    plt.close(fig)

    with open(f"{OUT}/report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report))
    print(f"\nFigures and report saved in '{OUT}'.")


if __name__ == "__main__":
    main()
