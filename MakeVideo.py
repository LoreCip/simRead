import os

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from .math_functions import find_zeros

from tqdm import tqdm
from p_tqdm import p_map

def determine_y_axis_range(data):
    min_value = np.min(data)
    max_value = np.max(data)

    if min_value > 0:
        y_min = 0
        y_max = max_value * 1.2  # Add some padding above the maximum value
    elif max_value < 0:
        y_min = min_value * 1.2  # Add some padding below the minimum value
        y_max = 0
    else:
        y_min = min_value * 1.2  # Add some padding below the minimum value
        y_max = max_value * 1.2  # Add some padding above the maximum value

    return [y_min, y_max]

def plotting(inputs):

    i, X, B, eps, rho, exp, t, m, path = inputs 

    DarkGray = '#233142'
    DarkBlue = '#455D7A'
    ClearRed = '#F95959'
    FullYellow = '#FACF5A'

    i, X, B, eps, rho, exp, t, m, path = inputs 

    matplotlib.rc('axes', edgecolor=DarkGray)
    fig, axes = plt.subplots(ncols=2, nrows=1, sharex = 'col', figsize=(10, 6))

    xlim = (0., 20)
    mask1 = (X > xlim[0]+0.1) 
    mask2 = (X < xlim[1])
    skipped = len(mask1) - sum(mask1)
    mask = mask1 & mask2

    ## PLOT 1
    # Epsilon
    axes[0].plot(X, eps, color = DarkBlue, zorder = 10)
    axes[0].set_xlim(xlim)
    axes[0].set_ylabel(r"$\epsilon$", rotation=0, color=DarkBlue)
    axes[0].set_ylim(determine_y_axis_range(eps[mask]))
    # B
    ax2 = axes[0].twinx() 
    ax2.set_ylabel('B', color = ClearRed, rotation=0)  
    ax2.get_yaxis().set_label_coords(1.13,0.5)
    ax2.plot(X, B, color = ClearRed)
    ax2.tick_params(axis='y')     
    ax2.set_ylim(determine_y_axis_range(B[mask]))

    ## PLOT 2
    # Rho
    axes[1].plot(X, rho, color = DarkBlue, zorder = 10)
    axes[1].set_xlim(xlim)
    axes[1].set_ylabel(r"$\rho$", rotation=0, color=DarkBlue)
    axes[1].set_ylim(determine_y_axis_range(rho[mask]))
   
    # Expansion
    ax2 = axes[1].twinx()
    ax2.set_ylabel(r'$\theta_+$', color = ClearRed, rotation=0)  
    ax2.get_yaxis().set_label_coords(1.12,0.5)
    ax2.plot(X, exp, color = ClearRed, zorder = 0)
    ax2.tick_params(axis='y')     
    ax2.set_ylim(determine_y_axis_range(exp[mask]))

    idxs = find_zeros(X, exp)
    if len(idxs) != 0: 
        plt.scatter(idxs[:,0], idxs[:, 1], color = DarkGray, marker = '.')

    axes[1].set_zorder(ax2.get_zorder()+1)
    axes[1].set_frame_on(False)

    textstr = "Time: " + str(np.round(t , 3))
    props = dict(boxstyle='round', facecolor=FullYellow, alpha=0.5)
    axes[1].text(0.8, 0.95, textstr, transform=axes[1].transAxes, fontsize=14, verticalalignment='top', horizontalalignment='center', bbox=props)

    for ii in range(2):
        axes[ii].set_xlabel('R')
    
    plt.tight_layout()

    plt.savefig(os.path.join(path, f'frame{i:07d}.png'))
    plt.close()

def ProduceInputs(sim, last_line, pline, path):

    it = np.array(range(last_line, last_line + pline))
    it = it[it < sim.niter]
    if len(it) == 0: return None

    X = sim.xgrid
    Bs = sim.get(it, 'B')
    ebs = sim.get(it, 'e^b')
    rhos = sim.get(it, 'rho')
    ts = sim.get(it, 't')
    exp = sim.get(it, 'expansion')
    m = sim.mass

    return [(last_line + i, X, Bs[i], ebs[i], rhos[i], exp[i], ts[i], m, str(path)) for i in range(0, len(it))]

def GenerateVideo(sim, ppath, n_partitions = 1, fps=30, verbose='quiet'):

    n_partitions = int(n_partitions)
    if n_partitions <= 0:
        n_partitions = 1

    Path(f"{os.path.join(ppath,'.frames')}").mkdir(parents=True, exist_ok=True)
    
    last_line = 0
    p_line = int(np.ceil( sim.niter / n_partitions))
    for j in range(n_partitions):
        inputs = ProduceInputs(sim, last_line, p_line, os.path.join(ppath,'.frames'))
        if inputs == None: break

        try:
            p_map(plotting, inputs)
        except:
            print('Parallel frames creation failed. Switching to serial.')
            for i, inp in tqdm(enumerate(inputs)):
                plotting(inp)

        last_line += p_line

    os.system(f"ffmpeg -framerate {fps} -loglevel {verbose} -pattern_type glob -i '{os.path.join(ppath,'.frames/frame*.png')}' -c:v libx264 -pix_fmt yuv420p '{os.path.join(ppath,'out_video.mp4')}'")
    os.system(f"rm -rf {os.path.join(ppath,'.frames')}")

    
