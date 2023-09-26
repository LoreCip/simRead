import os
import re
from pathlib import Path

import warnings
import errno

import h5py

import numpy as np
import scipy.optimize as op

import matplotlib.pyplot as plt

from . import MakeVideo
from . import math_functions as mf

class _Sims():

    def __init__(self, _path):

        self.check(_path)

        self.root = _path
        self.sims = {}

        for path in Path(self.root).rglob('*.h5'):

            match = re.search(r"id[\d]_m[\d.]+_dx[\d.]+_xMax[\d.]+_tf[\d.]+_r[\d.]+_a0[\d.]+(?:_corrected(?:_finer)?)?", str(path))
            if match:
                key = match.group()

            self.sims[key] = _Sim(str(path)[:-18], _name=key)

        self.keylist = list(self.sims.keys())


    def check(self, _path):
        if not os.path.isdir(_path):
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), _path)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.sims[key]
        elif isinstance(key, int):
            return self.sims[self.keylist[key]]
        elif isinstance(key, float):
            key = int(key)
            return self.sims[self.keylist[key]]
        elif isinstance(key, list) or isinstance(key, np.ndarray):
            return [self.__getitem__(k) for k in list(key)]
        elif isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or len(self.keylist)
            step = key.step or 1
            return [self.__getitem__(k) for k in range(start, stop, step)]

class _Sim():

    def __init__(self, _path, _name=None):

        self.check(_path)

        self.path = _path
        self.outpath = os.path.join(_path, 'outputs/output.h5')

        if _name == None:
            match = re.search(r"id[\d]_m[\d.]+_dx[\d.]+_xMax[\d.]+_tf[\d.]+_r[\d.]+_a0[\d.]+(?:_corrected(?:_finer)?)?", str(_path))
            if match:
                self.name = match.group()
        else:
            self.name = _name

        data = self.readParameters()
        self.id      = int(data[0])
        self.simtime = data[1]
        self.r0      = data[2]
        self.a0      = data[3]
        self.mass    = data[4]
        self.order   = int(data[5])
        self.xmax    = data[6]
        self.dx      = data[7]
        self.dSave   = int(data[8])
        self.dOut    = int(data[9])

        self.hor_loc = 2 * self.mass
        self.xb = self.xbounce()

        self.outfile = h5py.File(self.outpath, 'r')
        self.iterations = self.sort_groups()
        self.niter = len(self.iterations)

        self.xgrid = self.__getitem__('Xgrid')

        self.valid_keys = ['t', 'rho', 'e^b', 'B', 'beta', 'dt', 'theta', 'expansion']
        
    def __getitem__(self, key):
        item = self.outfile[key]
        if isinstance(item, h5py.Group):
            return GroupWrapper(item, self.xgrid)
        else:
            return item[()]

    def __str__(self):
        return f'Simulation: {self.name}'

    def check(self, _path):
        if not os.path.isdir(_path):
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), _path)

        if not os.path.isfile(os.path.join(_path, 'ParameterFile.par')):
            raise FileNotFoundError(os.path.join(_path, 'ParameterFile.par'))

        if not os.path.isfile(os.path.join(_path, 'outputs/output.h5')):
            raise FileNotFoundError(os.path.join(_path, 'outputs/output.h5'))

    def readParameters(self):
        comment_char = "/*"
        # Load the file using genfromtxt, ignoring comments
        return np.genfromtxt(
                                os.path.join(self.path, 'ParameterFile.par'),
                                comments=comment_char,
                                dtype = float
                            )

    def sort_groups(self):
        group_names = list(self.outfile.keys())
        rl = []
        for i in range(len(group_names)):
            try:
                int(group_names[i])  # Try converting the element to an integer
            except ValueError:
                rl.append(i)  # Remove the element if it cannot be converted
        for r in rl[::-1]:
            group_names.pop(r)

        return sorted(group_names, key=int)

    def get(self, iteration, item):

        if item not in self.valid_keys:
            raise IndexError("Invalid key: {} - Accepted keys: {}".format(item, self.valid_keys))

        if isinstance(iteration, int):
            indices = [iteration]
        elif isinstance(iteration, list) or isinstance(iteration, np.ndarray):
            indices = iteration
        elif isinstance(iteration, slice):
            start = iteration.start or 0
            stop = iteration.stop or self.niter
            step = iteration.step or 1
            indices = range(start, stop, step)
        else:
            try:
                iteration = int(iteration)
            except ValueError:
                raise ValueError("Invalid iteration value: {} - It should be an integer or convertible to one.")
            indices = [iteration]

        result = [self.__getitem__(self.iterations[i])[item] if not item in ['theta', 'expansion'] \
                                        else self.compute_expansion(iteration=i) for i in indices]

        if len(result) == 1:
            return result[0]
        return result

    def get_at_time(self, time, item, find_closest = True):
        
        if not isinstance(time, float):
            time = float(time)
        if time > self.simtime:
            raise ValueError("Requested time: {} - Maximum time available {}".format(time, self.simtime))
        if item not in self.valid_keys:
            raise IndexError("Invalid key: {} - Accepted keys: {}".format(item, self.valid_keys))

        # Max iteration possible if all timesteps where done with the largest timestep possible dt = dx/100
        max_iter = str(int(100 * time / self['dx']))
        # Find its index and add some iterations for safety
        try:
            idx = self.iterations.index(max_iter) + 20
        except ValueError as e:
            if not find_closest:
                raise e
            else:
                closest_idx = str(min(self.iterations, key=lambda x: abs(int(x) - int(max_iter))))
                idx = self.iterations.index(str(int(closest_idx))) + 20

        mf.find_i(idx, time, self)

        if not item in ['theta', 'expansion']:
            return self.__getitem__(self.iterations[i])[item]
        else:
            return self.compute_expansion(iteration = i )

    def xbounce(self):
        
        try:
            return op.brentq(mf.fxb, 
                        1e-6, 3*(self.hor_loc)**(1/3), 
                        args = (self.r0**2 / self.a0**2, self.hor_loc))
        except ValueError:
            psamp = np.linspace(self.xgrid[0], self.xgrid[-1], 1000)
            sampl = mf.fxb(psamp, self.r0**2 / self.a0**2, self.hor_loc)
            for i in range(1000):
                if sampl[i] * samp[i+1] <= 0:
                    break
            return op.brentq(mf.fxb, 
                        self.xgrid[i], self.xgrid[i+1], 
                        args = (self.r0**2 / self.a0**2, self.hor_loc))

    def find_timeout(self):

        BH_present = False
        never = True
        for it in range(self.niter):

            horizons = mf.find_zeros(self.xgrid, self.compute_expansion(iteration=it))
            if not BH_present and len(horizons) != 0:
                never = False
                BH_present = True
                t_formation = self.get(it, 't')
            elif BH_present and len(horizons) == 0:
                BH_present = False
                t_exit = self.get(it, 't')
                break

        if not never:
            return t_exit - t_formation
        else:
            return np.NaN

    def compute_expansion(self, iteration=None, time=None):
        """
                theta(i) = 1_RK - x(i)**2 * sin(2_RK * B(i) / x(i)**2)**2 / ( 4_RK * (1_RK + E(i)) )
        """

        if iteration == None and time != None:
            B = self.get_at_time(time, 'B')
            E = self.get_at_time(time, 'e^b')
        elif time == None and iteration != None:    
            B = self.get(iteration, 'B')
            E = self.get(iteration, 'e^b')
        elif iteration == None and time == None:
            raise ValueError("Both iteration and time are missing. Please provide either iteration or time.")
        else:
            warnings.warn("Both iteration and time given. Using iteration.")
            B = self.get(iteration, 'B')
            E = self.get(iteration, 'e^b')

        return mf.expan(self.xgrid, B, E)

    def make_video(self, batches=0, video_path=None, fps=30, verbose='quiet'):

        if video_path == None:
            video_path = os.path.join(self.path, 'video')
        Path(video_path).mkdir(parents=True, exist_ok=True)

        if batches == 0:
            batches = int(np.ceil(self.niter / 5000))

        if not isinstance(fps, int): fps = int(fps)
        if fps < 0: fps = 30

        if isinstance(verbose, int) and verbose%8 != 0:
            raise ValueError("Verbose should be -8 to not print any output, 0 to only show fatal errors or a multiple of 8 up to 56. Each level gives more details.")
        allowed = ["quiet", "panic", "fatal", "error", "warning", "info", "verbose", "debug", "trace"]
        if isinstance(verbose, str) and not verbose in allowed:
            raise ValueError(f"Verbose should be one of {allowed}.")

        MakeVideo.GenerateVideo(self, video_path, batches, fps, verbose)

    def plot(self, Y, X=None, iteration=None, time=None, color=None, linestyle=None, xrange=None, yrange=None, xlabel=None, ylabel=None, title=None, printTime=True, showHor=False, return_handles=False, savefig=False, name=None, savepath='.'):

        if iteration == None and time != None:
            qtt = self.get_at_time(time, Y)
            t = time
        elif time == None and iteration != None:    
            qtt = self.get(iteration, Y)
            t = self.get(iteration, 't')
        elif iteration == None and time == None:
            raise ValueError("Both iteration and time are missing. Please provide either iteration or time.")
        else:
            warnings.warn("Both iteration and time given. Using iteration.")
            qtt = self.get(iteration, Y)
            t = self.get(iteration, 't')

        if X == None:
            xgrid = self.xgrid
        else:
            xgrid = X

            if len(xgrid) != len(self.__getitem__('Xgrid')):
                raise ValueError("Length of provided grid is not equal to the length of the {} array.".format(item))

        if color != None:
            c = color
        else:
            c = 'k'
        if linestyle != None:
            ls = linestyle
        else:
            ls = '-'

        fig, ax = plt.subplots()
        plt.plot(xgrid, qtt, color = c, linestyle = ls)

        if showHor:
            plt.axvline(self.hor_loc, color = 'k', linestyle='--', alpha = 0.5)

        if xrange != None:
            if len(xrange) != 2:
                raise ValueError("Parameter xrange must be a tuple or list of lenght 2.")
            plt.xlim(xrange[0], xrange[1])
            
            if yrange == None:
                mask = (xgrid >= xrange[0]) & (xgrid <= xrange[1])
                bot = np.min(qtt[mask])
                top = np.max(qtt[mask])

                if np.abs(bot) > np.abs(top):
                    plt.ylim(bot*1.1, 0)
                else:
                    plt.ylim(0, top*1.1)

        if yrange != None:
            if len(yrange) != 2:
                raise ValueError("Parameter yrange must be a tuple or list of lenght 2.")
            plt.ylim(yrange[0], yrange[1])

        if xlabel != None:
            xl = str(xlabel)
        else:
            xl = 'Position'
        plt.xlabel(xl)
        
        if ylabel != None:
            yl = str(ylabel)
        else:
            yl = Y
        plt.ylabel(yl)

        if title != None:
            tl = str(title)
        else:
            tl = str(Y)
        plt.title(tl)

        if printTime:
            textstr = "Time: " + str(np.round(t, 3))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            plt.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

        if savefig:
            if name == None:
                name = f'{str(Y)}_time{np.round(t,3)}.png'
            if not name[:-4] in ['.png', '.jpg']:
                name = name + '.png'
            savepath = os.path.join(savepath, name)

            plt.savefig(savepath, format='png')

        if return_handles:
            return fig, ax

class GroupWrapper:
    def __init__(self, group, xgrid):
        self.group = group
        self.conv = False
        self.xgrid = xgrid

    def __getitem__(self, key):

        if key == 'beta':
            key = 'B'
            self.conv = True

        if not self.conv:
            return self.group[key][()]
        else:
            return self.group[key][()] / self.xgrid**2