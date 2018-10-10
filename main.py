from sys import argv
from os import listdir
from numpy import *
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

def get_data(filename, dir):
    try:
        r = get_r(filename)
        t = get_theta(filename)
        e = get_energy(dir + "/" + filename)
        return (r,t,e)
    except:
        print("Failed to get data point for " + filename)

def fit_data(data):
    r_eq, t_eq, e_eq = min(data, key = lambda tup: tup[2])
    print(r_eq)
    print(t_eq)
    print(e_eq)

def plot(data):
    #Format data so it can be used by matplotlib
    es = { (r,t): e for (r,t,e) in data }

    rs, ts, _ = zip(*data)    
    rs = sorted(rs)
    ts = sorted(ts)

    R,T = meshgrid(rs,ts)
    E = np.array([ np.array([es[(r,t)] for r in rs]) for t in ts])

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.contour(R,T,E)
    plt.show()

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

    data = [get_data(filename, dir) for filename in file_list if get_data(filename, dir) != None]

    fit_data(data)
    plot(data)

main()
