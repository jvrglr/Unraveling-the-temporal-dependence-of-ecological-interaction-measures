# %%
%reset
import matplotlib
from scipy.optimize import curve_fit
import numpy as np
from numpy import exp as exp
from numpy import log as log
from numpy import sqrt as sqrt
import networkx as nx
import matplotlib.pyplot as plt
import time
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.patches as patches
matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size':12})
matplotlib.rcParams['axes.linewidth']=1 #Grosor del marco (doble del standard)
def Plot_bonito(xlabel=r" $ x$",ylabel=r"$ y$",label_font_size=15,ticks_size=12,y_size=2.4,x_size=3.2):
    plt.figure(figsize=(x_size,y_size))
    plt.tick_params(labelsize=24)
    plt.xlabel(xlabel,fontsize=label_font_size)
    plt.ylabel(ylabel,fontsize=label_font_size)
    plt.xticks(fontsize=ticks_size)
    plt.yticks(fontsize=ticks_size)
    plt.locator_params(axis="both", nbins=5,tight=True)
def axis_bonito(ax,xlabel=r" $ x$",ylabel=r"$ y$",label_font_size=12,ticks_size=10):
    ax.set_xlabel(xlabel,fontsize=label_font_size)
    ax.set_ylabel(ylabel,fontsize=label_font_size)
    ax.tick_params(axis="x", labelsize=ticks_size)
    ax.tick_params(axis="y", labelsize=ticks_size)
    # For two column figure that fits in one column of a two-columns paper:
        #fig, (ax1,ax2) = plt.subplots(1, 2,figsize=(4.2,1.8))
    # For three rows figure without space between plots:
        #fig, (ax1,ax2,ax3) = plt.subplots(3, 1,sharex=True,figsize=(2,3.6))
        #fig.tight_layout()
       # #-----plots--------
        #plt.subplots_adjust(hspace=0)
#%% Functions

# Compute relative gradients
def compute_gradient(data):
    x = data["x"].to_numpy()
    t = data["t"].to_numpy()
    dx = np.diff(x)
    dt = np.diff(t)
    if np.any(dt < 0):
        print("Warning: Time values are not strictly increasing.")
    return dx / dt / x[:-1]
# %%
Plot_bonito(xlabel=r"$t\,\,(\mathrm{days}) $",ylabel=r"$C\,\,(\mathrm{m}^{-2})$")
names = "M S P".split(" ")
markers = "o s ^".split(" ")
colors = "dodgerblue salmon mediumseagreen".split(" ")
exp_decay = [-0.15,-0.3,-0.2]
for ii,name in enumerate(names):
    d = pd.read_csv("data/Grasshoppers/"+name+"_alone.csv",names=["t","x"],dtype=float,skiprows=1)
    ts = d["t"].to_numpy()
    xs = d["x"].to_numpy()/0.16
    print("ii=",ii,"name=",name,"C(0)=",xs[0])
    
    #Data
    plt.scatter(ts,xs,s=8,label=name,marker=markers[ii],color=colors[ii])
    plt.plot(ts,xs,color=colors[ii],ls="--",alpha=0.2)

    #Exponential fit
    c = xs[-1]
    a = xs[0] - c
    plt.plot(ts,a*exp(exp_decay[ii]*ts)+c,color=colors[ii],lw=1)

    # Biological fit
    data = np.loadtxt("data/Grasshoppers/traj_"+names[ii]+"_alone.dat").T
    ts = data[0]
    Cs = data[1:4]
    Rs = data[4:]
    dt = ts[1]-ts[0]
    plt.plot(ts,Cs[ii],ls="-",color=colors[ii],alpha=0.3,lw=5)

plt.yscale("log")
plt.xticks([0,10,20,30])
plt.show();plt.close()

# %% Computing interactions over data
Plot_bonito(xlabel=r"$t\,\,(\mathrm{days}) $",ylabel=r"$\mathcal{M} \,\,(\mathrm{m}^{-2}\mathrm{day}^{-1})$",x_size=4,y_size=2.4)
baseline_names = ["M", "S", "P"]
perturbed_names = ["MF", "SC", "PN"]
colors = ["dodgerblue", "salmon", "mediumseagreen", "black", "mediumpurple", "darkorange"]
alphas = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
labels = [r"$\mathrm{S} \to \mathrm{M}$",r"$\mathrm{P} \to \mathrm{M}$",r"$\mathrm{M} \to \mathrm{S}$",r"$\mathrm{P} \to \mathrm{S}$",r"$\mathrm{M} \to \mathrm{P}$",r"$\mathrm{S} \to \mathrm{P}$"]
markers = ["o", "s", "D", "^", "v", "P"]  # circle, square, diamond, triangle_up, triangle_down, plus-filled

count = -1
for baseline_index,baseline_name in enumerate(baseline_names):
    for perturbed_index,perturbed_name in enumerate(perturbed_names):
        if baseline_index == perturbed_index:
            continue
        else:
            count += 1

        # Load data
        baseline_path = f"data/Grasshoppers/{baseline_name}_alone.csv"
        perturbed_path = f"data/Grasshoppers/{baseline_name}_w_{perturbed_name}.csv"

        col_names = ["t", "x"]
        baseline_data = pd.read_csv(baseline_path, names=col_names, dtype=float, skiprows=1).sort_values("t")
        baseline_data["t"] = baseline_data["t"].round().astype(int)

        perturbed_data = pd.read_csv(perturbed_path, names=col_names, dtype=float, skiprows=1).sort_values("t")
        perturbed_data["t"] = perturbed_data["t"].round().astype(int)
        if len(perturbed_data) != len(baseline_data):
            print("Warning: Data lengths do not match.")

            # Round times for approximate matching
            baseline_times = baseline_data["t"].round().astype(int)
            perturbed_data["t_rounded"] = perturbed_data["t"].round().astype(int)

            if len(perturbed_data) > len(baseline_data):
                # Keep only rows in perturbed_data where rounded t is in baseline rounded t
                perturbed_data = perturbed_data[perturbed_data["t_rounded"].isin(baseline_times)].reset_index(drop=True)
            else:
                # Keep only rows in baseline_data where rounded t is in perturbed rounded t
                baseline_data["t_rounded"] = baseline_data["t"].round().astype(int)
                perturbed_times = perturbed_data["t_rounded"]
                baseline_data = baseline_data[baseline_data["t_rounded"].isin(perturbed_times)].reset_index(drop=True)

            # Drop helper column after filtering
            perturbed_data.drop(columns=["t_rounded"], inplace=True, errors='ignore')
            baseline_data.drop(columns=["t_rounded"], inplace=True, errors='ignore')

                

        gradient_baseline = compute_gradient(baseline_data)
        gradient_perturbed = compute_gradient(perturbed_data)
        gradient_diff = gradient_perturbed - gradient_baseline

        # Plot
        time = baseline_data["t"].to_numpy()[:-1]
        color = colors[count]
        s = 20
        lw = 0.5
        plt.scatter(time, gradient_diff, color=color,marker=markers[count],s=s,label=labels[count],alpha=alphas[count],edgecolors="none")
        plt.plot(time, gradient_diff, color=color, linestyle="--", alpha=0.5,lw = lw)
        plt.plot(time, np.zeros_like(time), color="black", linestyle="-", alpha=0.2, linewidth=0.5)
plt.legend(loc="upper right", fontsize=10, frameon=False)
plt.legend().remove()  # remove legend from the current plot
plt.savefig("figures/Grasshoppers/gradient_diff.pdf", bbox_inches='tight',transparent=True)

# Save legend
legend_handles = [
    plt.Line2D([0], [0], marker=markers[i], color='w', label=labels[i],
               markerfacecolor=colors[i], markersize=6, linestyle='None')
    for i in range(len(labels))
]

# Create legend-only figure
fig_legend = plt.figure(figsize=(6, 1))
ax = fig_legend.add_subplot(111)
ax.axis('off')  # no axes

legend = ax.legend(
    handles=legend_handles,
    loc='center',
    ncol=3,
    frameon=False,
    fontsize=10
)

# Save legend-only figure
fig_legend.savefig("figures/Grasshoppers/legend_only.pdf", bbox_inches='tight')


# %% Computing interactions over data (with subplots + bootstrap)

# ---------------- Figure setup ----------------
fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey=True)
axes = axes.flatten()

baseline_names = ["M", "S", "P"]
perturbed_names = ["MF", "SC", "PN"]

colors = ["dodgerblue", "salmon", "mediumseagreen",
          "mediumpurple", "darkorange", "black"]

labels = [r"$\mathrm{S} \to \mathrm{M}$",
          r"$\mathrm{P} \to \mathrm{M}$",
          r"$\mathrm{M} \to \mathrm{S}$",
          r"$\mathrm{P} \to \mathrm{S}$",
          r"$\mathrm{M} \to \mathrm{P}$",
          r"$\mathrm{S} \to \mathrm{P}$"]

markers = ["o", "s", "D", "^", "v", "P"]

n_boot = 1000
count = -1

# ==========================================================
# Main loop
# ==========================================================

for baseline_index, baseline_name in enumerate(baseline_names):
    for perturbed_index, perturbed_name in enumerate(perturbed_names):

        if baseline_index == perturbed_index:
            continue

        count += 1
        ax = axes[count]

        # ---------------- Load data ----------------
        baseline_path = f"data/Grasshoppers/{baseline_name}_alone.csv"
        perturbed_path = f"data/Grasshoppers/{baseline_name}_w_{perturbed_name}.csv"

        col_names = ["t", "x"]

        baseline_data = pd.read_csv(
            baseline_path, names=col_names, dtype=float, skiprows=1
        ).sort_values("t")

        perturbed_data = pd.read_csv(
            perturbed_path, names=col_names, dtype=float, skiprows=1
        ).sort_values("t")

        # ---------------- Align time points ----------------
        common_times = np.intersect1d(
            baseline_data["t"].round().astype(int),
            perturbed_data["t"].round().astype(int)
        )

        baseline_data = baseline_data[
            baseline_data["t"].round().astype(int).isin(common_times)
        ].reset_index(drop=True)

        perturbed_data = perturbed_data[
            perturbed_data["t"].round().astype(int).isin(common_times)
        ].reset_index(drop=True)

        # ---------------- Compute gradients ----------------
        gradient_baseline = compute_gradient(baseline_data)
        gradient_perturbed = compute_gradient(perturbed_data)

        # Ensure equal length
        min_len = min(len(gradient_baseline), len(gradient_perturbed))
        gradient_baseline = gradient_baseline[:min_len]
        gradient_perturbed = gradient_perturbed[:min_len]

        gradient_diff = gradient_perturbed - gradient_baseline

        time = baseline_data["t"].to_numpy()[:min_len]
        time = time[:-1]  # gradients reduce one time point

        gradient_diff = gradient_diff[:len(time)]

        # =====================================================
        # Residual bootstrap
        # =====================================================

        # Simple smoothing (moving average)
        window = 3
        smooth = pd.Series(gradient_diff).rolling(
            window, center=True, min_periods=1
        ).mean().to_numpy()

        residuals = gradient_diff - smooth

        boot_samples = np.zeros((n_boot, len(gradient_diff)))

        for b in range(n_boot):
            resampled_res = np.random.choice(
                residuals, size=len(residuals), replace=True
            )
            boot_samples[b] = smooth + resampled_res

        lower = np.percentile(boot_samples, 2.5, axis=0)
        upper = np.percentile(boot_samples, 97.5, axis=0)

        # =====================================================
        # Plot
        # =====================================================

        color = colors[count]

        # Confidence band
        ax.fill_between(time, lower, upper,
                        color=color, alpha=0.2, linewidth=0)

        # Interaction points
        ax.scatter(time, gradient_diff,
                   color=color,
                   marker=markers[count],
                   s=20,
                   edgecolors="none")

        # Line
        ax.plot(time, gradient_diff,
                color=color,
                linestyle="--",
                alpha=0.6)

        # Zero line
        ax.axhline(0, color="black",
                   linestyle="-",
                   linewidth=0.5,
                   alpha=0.3)

        ax.set_title(labels[count], fontsize=10)

# ==========================================================
# Global labels
# ==========================================================

fig.supxlabel(r"$t\,\,(\mathrm{days})$")
fig.supylabel(r"$\mathcal{M}\,\,(\mathrm{m}^{-2}\mathrm{day}^{-1})$")

plt.tight_layout()

plt.savefig("figures/Grasshoppers/gradient_diff_subplots.pdf",
            bbox_inches='tight', transparent=True)
# %%
