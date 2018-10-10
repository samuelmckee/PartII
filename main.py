from sys import argv
from os import listdir
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def get_theta(filename):
    theta = filename.split("theta")[1][:-4]
    return float(theta)


def get_r(filename):
    r = filename.split("theta")[0].split("r")[-1]
    return float(r)


def get_energy(file_path):
    f = open(file_path, "r")
    for line in f:
        if "SCF Done" in line:
            l = line.split()
            return float(l[4])

def main():
    try:
        dir = argv[1]
    except IndexError:
        print("Please enter directory as command line containing output files as argument.")
        return #Exit program

    file_list = [str for str in listdir(dir) if str.endswith(".out")]

    if len(file_list) == 0:
        print("No .out files found in given directory")
        return #Exit program

    rs = [get_r(filename) for filename in file_list]
    ts = [get_theta(filename) for filename in file_list]
    es = [get_energy(dir + "/" + filename) for filename in file_list]

    es = { (rs[i],ts[i]): es[i] for i in range(len(rs)) }

    rs = sorted(set(rs))
    ts = sorted(set(ts))

    R,T = np.meshgrid(rs,ts)
    E = np.array([ np.array([es[(r,t)] for r in rs]) for t in ts])

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(R,T,E)
    plt.show()


main()
