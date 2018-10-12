from sys import argv
from os import listdir
from datetime import datetime
from scipy.interpolate import griddata
from numpy import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


def main():
    try:
        dir = argv[1]
    except IndexError:
        print('Please enter directory containing output files as command '
            'line argument.')
        return #Exit program

    file_list = [str for str in listdir(dir) if str.endswith('.out')]

    if len(file_list) == 0:
        print('No .out files found in given directory')
        return #Exit program

    data = [get_data(filename, dir) for filename
        in file_list if get_data(filename, dir) != None]

    print('Loaded ' + str(len(data)) + ' files')

    dt = datetime.now()
    out = '%d-%02d-%02d_%02d:%02d:%02d' % (dt.year,
                                 dt.month,
                                 dt.day,
                                 dt.hour,
                                 dt.minute,
                                 dt.second)
    plot(data, 'FIGURE_' + out + '.png')


# Functions for importing data from files
def get_theta(filename):
    theta = filename.split('theta')[1][:-4]
    return float(theta)


def get_r(filename):
    r = filename.split('theta')[0].split('r')[-1]
    return float(r)


def get_energy(file_path):
    f = open(file_path, 'r')
    for line in f:
        if 'SCF Done' in line:
            l = line.split()
            return float(l[4])


def get_data(filename, dir):
    try:
        r = get_r(filename)
        t = get_theta(filename)
        e = get_energy(dir + '/' + filename)
        return (r,t,e)
    except:
        print('Failed to get data point for ' + filename)


# Functions for processing data
def plot(data,out):
    # Prepare points to be plotted
    # Values are interpolated to give smoother colour gradient on plot
    r,t,e = zip(*data)
    r_space = linspace(min(r),max(r),100)
    t_space = linspace(min(t),max(t),100)
    grid_r, grid_t = meshgrid(r_space, t_space)
    grid_e = griddata((r, t), e, (grid_r, grid_t), method='cubic')

    #Plot points and save figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(30,285)
    ax.set_xlabel('r / Angstroms')
    ax.set_ylabel('theta / degrees')
    ax.set_zlabel('Energy / Hartrees')
    ax.grid(False)
    ax.plot_surface(grid_r,grid_t,grid_e, rstride=1, cstride=1,
        cmap='inferno')
    plt.savefig(out, dpi = 200)
    print('Saved surface plot as ' + out)

main()
