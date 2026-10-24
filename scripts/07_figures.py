"""Figures for the README (outputs/figures)."""
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
TAB = ROOT / "outputs" / "tables"
GIS = ROOT / "outputs" / "gis"
FIG = ROOT / "outputs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
AQ_ORDER = ["Ogallala (High Plains)", "Edwards-Trinity (Plateau)", "Carrizo-Wilcox", "Gulf Coast"]
AQ_COL = dict(zip(AQ_ORDER, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
COMP_COL = {"TWSA": "#2a78d6", "Soil moisture": "#eb6834", "Reservoirs": "#1baf7a", "GWSA": "#4a3aa7"}
DIVERGE = LinearSegmentedColormap.from_list("rdgybu", ["#a8302f", "#e34948", "#f0efec", "#3987e5", "#1c5cab"])
DROUGHTS = [("2010-10", "2011-12"), ("2022-03", "2023-09")]
GAP = ("2017-07", "2018-05")

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10,
    "lines.linewidth": 1.4, "legend.frameon": False})


def dt(m):
    return pd.PeriodIndex(m, freq="M").to_timestamp() + pd.Timedelta(days=14)


def shade(ax):
    for a, b in DROUGHTS:
        ax.axvspan(pd.Timestamp(a), pd.Timestamp(b) + pd.offsets.MonthEnd(), color="#f2d3c9", alpha=.55, lw=0, zorder=0)
    ax.axvspan(pd.Timestamp(GAP[0]), pd.Timestamp(GAP[1]) + pd.offsets.MonthEnd(), color="#e4e3df", alpha=.8, lw=0, zorder=0)


def fig_map():
    m = gpd.read_file(GIS / "grace_mascon_trends_texas.gpkg")
    aq = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="aquifers")
    tx = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="texas")
    lim = np.ceil(np.nanpercentile(np.abs(m.trend_mm_yr), 99))
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    m.plot(column="trend_mm_yr", ax=ax, cmap=DIVERGE, norm=TwoSlopeNorm(0, -lim, lim),
           edgecolor="white", linewidth=0.3, legend=True,
           legend_kwds={"label": "TWSA trend, Apr 2002 to Aug 2026 (mm/yr)", "shrink": .7})
    sig = m[m.p_value < 0.05]
    ax.scatter(sig.lon_c, sig.lat_c, s=4, color=INK, zorder=4)
    tx.boundary.plot(ax=ax, color=INK, linewidth=0.9)
    for name, g in aq.groupby("aquifer"):
        g.boundary.plot(ax=ax, color=INK2, linewidth=0.7, linestyle="--")
        c = g.geometry.iloc[0].representative_point()
        ax.annotate(name, (c.x, c.y), fontsize=7.5, ha="center", color=INK,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=.75))
    ax.set_xlim(-107, -93.2)
    ax.set_ylim(25.5, 36.8)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("GRACE/GRACE-FO TWSA trend per CSR RL06.3 mascon\n"
                 "(dots: p < 0.05; dashed: TWDB major aquifers)", loc="left")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_mascon_twsa_trend_map.png", dpi=170)
    plt.close(fig)


def fig_components():
    df = pd.read_csv(PROC / "storage_components_aquifer_monthly.csv")
    fig, axes = plt.subplots(4, 1, figsize=(8.5, 10), sharex=True)
    for ax, aq in zip(axes, AQ_ORDER):
        d = df[df.aquifer == aq].copy()
        # connect across isolated missing months; break the line only at the
        # 11-month GRACE to GRACE-FO gap
        d = d.set_index("month")
        d.loc["2017-12"] = np.nan
        d = d.sort_index()
        x = dt(d.index)
        shade(ax)
        for col, lab in [("twsa_mm", "TWSA"), ("sma_mm", "Soil moisture"), ("rsa_mm_in", "Reservoirs"), ("gwsa_mm", "GWSA")]:
            ax.plot(x, d[col], color=COMP_COL[lab], lw=2.0 if lab in ("TWSA", "GWSA") else 1.2, label=lab)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(aq, loc="left")
        ax.set_ylabel("Anomaly (mm)")
    axes[0].legend(ncol=4, loc="lower left", fontsize=8)
    axes[-1].set_xlabel("Shaded: 2010-11 and 2022-23 droughts (pink), GRACE to GRACE-FO gap (gray)")
    fig.suptitle("Storage anomalies relative to 2004 to 2009, GWSA = TWSA - soil moisture - reservoirs", x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_storage_components.png", dpi=160)
    plt.close(fig)


def fig_stl():
    s = pd.read_csv(PROC / "stl_components_monthly.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [2.4, 1]})
    for aq in AQ_ORDER:
        d = s[(s.aquifer == aq) & (s.component == "GWSA")]
        x = dt(d.month)
        axes[0].plot(x, d.trend, color=AQ_COL[aq], lw=2, label=aq)
        clim = d.assign(m=pd.PeriodIndex(d.month, freq="M").month).groupby("m").seasonal.mean()
        axes[1].plot(clim.index, clim.values, color=AQ_COL[aq], lw=2, marker="o", ms=4, label=aq)
    shade(axes[0])
    axes[0].axhline(0, color=INK2, lw=.6)
    axes[0].set_title("STL trend of GWSA (period 12, outliers down-weighted)", loc="left")
    axes[0].set_ylabel("mm")
    axes[0].legend(fontsize=8, loc="lower left")
    axes[1].set_title("Mean seasonal component", loc="left")
    axes[1].set_xticks(range(1, 13))
    axes[1].set_xticklabels(list("JFMAMJJASOND"))
    axes[1].axhline(0, color=INK2, lw=.6)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_stl_gwsa.png", dpi=160)
    plt.close(fig)


def fig_drought():
    t = pd.read_csv(TAB / "table3_drought_response.csv", dtype={"drought": str})
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    for ax, comp in zip(axes, ["TWSA", "GWSA"]):
        d = t[t.component == comp]
        w = 0.38
        for i, ev in enumerate(["2011", "2022"]):
            e = d[d.drought == ev].set_index("aquifer").loc[AQ_ORDER]
            xs = np.arange(4) + (i - .5) * w
            ax.bar(xs, e.drawdown_km3, width=w - 0.03, color=["#86b6ef", "#1c5cab"][i], label=f"{ev} drought")
            for xx, v in zip(xs, e.drawdown_km3):
                ax.text(xx, v - 0.4, f"{v:.1f}", ha="center", va="top", fontsize=7, color=INK2)
        ax.set_xticks(range(4))
        ax.set_xticklabels([a.replace(" (", "\n(") for a in AQ_ORDER], fontsize=7.5)
        ax.axhline(0, color=INK2, lw=.6)
        ax.set_title(f"{comp}: drop from pre-drought year to drought minimum", loc="left", fontsize=9)
    axes[0].set_ylabel("Deseasonalised change (km3)")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_drought_drawdown.png", dpi=160)
    plt.close(fig)


def fig_wells():
    a = pd.read_csv(PROC / "annual_grace_vs_wells.csv").dropna(subset=["well_anom_m"])
    c = pd.read_csv(TAB / "table4_well_correlation.csv")
    fig, axes = plt.subplots(2, 4, figsize=(12, 6))
    for j, aq in enumerate(AQ_ORDER):
        d = a[a.aquifer == aq].sort_values("year")
        # top: z-scores on one axis, so no second y-scale is needed
        ax = axes[0, j]
        z = lambda v: (v - v.mean()) / v.std()  # noqa: E731
        ax.plot(d.year, z(d.gwsa_mm), color=COMP_COL["GWSA"], lw=2, marker="o", ms=3.5, label="GRACE GWSA")
        ax.plot(d.year, z(d.well_anom_m), color=INK2, lw=1.6, marker="s", ms=3.5, label="Well level (median)")
        ax.set_title(aq, loc="left", fontsize=9)
        if j == 0:
            ax.set_ylabel("Standardised anomaly (z)")
            ax.legend(fontsize=7.5, loc="lower left")
        ax = axes[1, j]
        ax.scatter(d.well_anom_m, d.gwsa_mm, s=22, color=COMP_COL["GWSA"], edgecolor="white", lw=.6, zorder=3)
        k = np.polyfit(d.well_anom_m, d.gwsa_mm, 1)
        xx = np.linspace(d.well_anom_m.min(), d.well_anom_m.max(), 10)
        ax.plot(xx, np.polyval(k, xx), color=INK2, lw=1)
        r = c[(c.aquifer == aq) & (c.grace_component == "GWSA") & (c.scale == "annual")].iloc[0]
        rt = c[(c.aquifer == aq) & (c.grace_component == "TWSA") & (c.scale == "annual")].iloc[0]
        ptxt = "p < 0.001" if r.p < 0.001 else f"p = {r.p:.2f}"
        ax.text(0.03, 0.97, f"GWSA r = {r.pearson_r:.2f} ({ptxt})\nTWSA r = {rt.pearson_r:.2f}\nn = {int(r.n)} years",
                transform=ax.transAxes, va="top", fontsize=7.5, zorder=5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=GRID, alpha=.9))
        ax.set_xlabel("Well-level anomaly (m)")
        if j == 0:
            ax.set_ylabel("GRACE GWSA, annual mean (mm)")
    fig.suptitle("GRACE groundwater storage anomaly vs TWDB GWDB well-level anomalies (annual, relative to 2004 to 2009)",
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIG / "fig5_grace_vs_wells.png", dpi=160)
    plt.close(fig)


def fig_trends():
    t = pd.read_csv(TAB / "table1_trends.csv")
    t = t[t.component.isin(["TWSA", "GWSA"])]
    periods = ["full 2002-04 to 2026-08", "GRACE 2002-04 to 2017-06", "GRACE-FO 2018-06 to 2026-08"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, p in zip(axes, periods):
        for k, comp in enumerate(["TWSA", "GWSA"]):
            d = t[(t.period == p) & (t.component == comp)].set_index("aquifer").loc[AQ_ORDER]
            ys = np.arange(4) + (k - .5) * 0.3
            ax.errorbar(d.trend_mm_yr, ys, xerr=[d.trend_mm_yr - d.ci95_lo, d.ci95_hi - d.trend_mm_yr],
                        fmt="o", ms=6, color=COMP_COL[comp], ecolor=COMP_COL[comp], elinewidth=1.4, capsize=0, label=comp)
        ax.axvline(0, color=INK2, lw=.7)
        ax.set_title(p, loc="left", fontsize=9)
        ax.set_xlabel("Trend (mm/yr), 95% CI (Newey-West)")
    axes[0].set_yticks(range(4))
    axes[0].set_yticklabels(AQ_ORDER)
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "fig6_trends_by_period.png", dpi=160)
    plt.close(fig)


def fig_leakage():
    t = pd.read_csv(TAB / "table5_leakage_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for aq in AQ_ORDER:
        d = t[t.aquifer == aq]
        ax.plot(d.mask_buffer_km, d.twsa_trend_mm_yr, color=AQ_COL[aq], marker="o", ms=6, lw=2, label=aq)
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks([-50, 0, 100, 200])
    ax.set_xlabel("Aquifer mask buffer (km; negative = interior only)")
    ax.set_ylabel("TWSA trend (mm/yr)")
    ax.set_title("Sensitivity of the TWSA trend to the averaging mask", loc="left")
    ax.legend(fontsize=7.5, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIG / "fig7_mask_sensitivity.png", dpi=160)
    plt.close(fig)


def fig_wellmap():
    w = pd.read_csv(PROC / "wells_annual_set.csv")
    aq = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="aquifers")
    tx = gpd.read_file(GIS / "aquifer_regions.gpkg", layer="texas")
    fig, ax = plt.subplots(figsize=(6.5, 6))
    tx.boundary.plot(ax=ax, color=INK, lw=.9)
    for name in AQ_ORDER:
        aq[aq.aquifer == name].plot(ax=ax, color=AQ_COL[name], alpha=.25, edgecolor="none")
        d = w[w.aquifer == name]
        ax.scatter(d.lon, d.lat, s=3, color=AQ_COL[name], label=f"{name} ({len(d)} wells)")
    ax.legend(fontsize=7.5, loc="lower left", markerscale=3)
    ax.set_title("TWDB GWDB wells in the annual comparison set (>= 15 years since 2002)", loc="left", fontsize=9)
    ax.grid(False)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    fig.tight_layout()
    fig.savefig(FIG / "fig8_well_locations.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig_map()
    fig_components()
    fig_stl()
    fig_drought()
    fig_wells()
    fig_trends()
    fig_leakage()
    fig_wellmap()
    print("figures written:", sorted(p.name for p in FIG.glob("*.png")))
