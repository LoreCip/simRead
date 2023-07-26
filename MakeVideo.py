import os

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from p_tqdm import p_map

########
## TO DO
########
# 1) check for existance of frame folder and change name
# 2) check for existance of out_video and change name
# 3) input name of out_video
# 4) option to keep frames

def determine_y_axis_range(data):
    min_value = np.min(data)
    max_value = np.max(data)

    if min_value > 0:
        y_min = 0
        y_max = max_value * 1.3  # Add some padding above the maximum value
    elif max_value < 0:
        y_min = min_value * 1.3  # Add some padding below the minimum value
        y_max = 0
    else:
        y_min = min_value * 1.3  # Add some padding below the minimum value
        y_max = max_value * 1.3  # Add some padding above the maximum value

    return [y_min, y_max]

def plotting(inputs):

        i, X, B, eps, rho, t, m, path = inputs 

        fig = plt.figure(figsize=(13, 10))
        fig.suptitle(f'Time: {np.round(t, 4)}', fontsize=14)
        ax1 = fig.add_subplot(3,1,1, adjustable='box')
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)
        ax1.plot(X, rho)
        ax2.plot(X, B)
        ax3.plot(X, eps)

        xlim=[0, 2*m*2.5]

        ax1.set(xlim=xlim, ylim=determine_y_axis_range(rho[100:-100]))
        ax2.set(xlim=xlim, ylim=determine_y_axis_range(B[100:-100]))
        ax3.set(xlim=xlim, ylim=determine_y_axis_range(eps[100:-100]))

        ax2.set_ylabel('B')
        ax3.set_ylabel(r'$\epsilon^b$')
        ax1.set_ylabel(r'$\rho$')

        ax3.set_xlabel('x')

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
    m = sim.mass

    return [(last_line + i, X, Bs[i], ebs[i], rhos[i], ts[i], m, str(path)) for i in range(0, len(it))]

def GenerateVideo(sim, ppath, n_partitions = 1):

    n_partitions = int(n_partitions)
    if n_partitions <= 0:
        n_partitions = 1

    Path(f"{os.path.join(ppath,'.frames')}").mkdir(parents=True, exist_ok=True)
    
    last_line = 0
    p_line = int(np.ceil( sim.niter / n_partitions))
    for j in range(n_partitions):
        inputs = ProduceInputs(sim, last_line, p_line, os.path.join(ppath,'.frames'))
        if inputs == None: break

        p_map(plotting, inputs)
        last_line += p_line

    os.system(f"ffmpeg -framerate 30 -loglevel quiet -pattern_type glob -i {os.path.join(ppath,'.frames/*.png')} -c:v libx264 -pix_fmt yuv420p {os.path.join(ppath,'video/out_video.mp4')}")
    os.system("rm -rf ./.frames")

    