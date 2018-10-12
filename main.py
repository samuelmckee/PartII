from sys import argv
from os import listdir
from math import *
from numpy import *

def main():
    try:
        dir = argv[1]
    except IndexError:
        print("Please enter directory as command line "
                    "containing output files as argument.")
        return # Exit program

    file_list = [str for str in listdir(dir) if str.endswith(".out")]

    if len(file_list) == 0:
        print("No .out files found in given directory")
        return # Exit program

    data = [get_data(filename, dir) for filename
                in file_list if get_data(filename, dir) != None]

    fit(data)


# Functions for importing data from files
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
        print("Failed to read data from " + filename)


# Functions for processing data

# Returns r^2 value for polynomial p which is fitted to xs and ys
def get_r2(p, xs, ys):
    y_fit = [polyval(p, x) for x in xs] # y values predicted by parameters
    y_avg = sum(ys)/len(ys)             # Average y value
    SStot = sum((ys - y_avg)**2)        # Sum of squares
    y_res = subtract(ys, y_fit)         # Residuals
    SSres = sum(y_res**2)               # Residual sum of squares

    return 1. - SSres / SStot

# Fits quadratic to xs and ys, reducing maximum displacement from equilibrium
# Until R^2 rises above minimum value
def fit_to_tol(xs, ys, x_eq, r2):
    p = polyfit(xs, ys, 2) # Fit data to quadratic
    range = max([abs(x-x_eq) for x in xs])
    step = range / 20. # Amount max displacement reduced by after each fit
    R2 = get_r2(p, xs, ys)

    while R2 < r2:     # R^2 below threshold
        #Reduce maximum displacement from equilibrium
        range -= step
        ys = [y for (i, y) in enumerate(ys) if abs(xs[i] - x_eq) < range]
        xs = [x for x in xs if abs(x - x_eq) < range]

        #Recalculate fit
        p = polyfit(xs, ys, 2)
        R2 = get_r2(p,xs,ys)

    print (len(xs))
    return p[0]

def fit(data):
    r_eq, t_eq, e_eq = min(data, key = lambda tup: tup[2])

    r_fit = [(r, t, e) for (r, t, e) in data if t == t_eq ]  # Keep t = t_eq
    r_fit.sort(key = lambda tup: tup[0])
    rs = [x[0] * 1.0e-10 for x in r_fit]      # Get r in m
    es = [x[2] * 4.35974e-18 for x in r_fit]  # Get E in J
    k_r = 2. * fit_to_tol(rs, es, r_eq * 1.0e-10, 0.95)

    t_fit = [(r, t, e) for (r, t, e) in data if r == r_eq]  # Keep r = r_eq
    t_fit.sort(key = lambda tup: tup[1])
    ts = [x[1] * pi / 180. for x in t_fit]    # Get theta in radians
    es = [x[2] * 4.35974e-18 for x in t_fit]  # Get E in J
    k_t = 2. * fit_to_tol(ts, es, t_eq * pi / 180., 0.95)

    mu_1 = 2.0 * 1.66e-27
    mu_2 = 0.5 * 1.66e-27

    nu_1 = sqrt(k_r / mu_1) / (2. * pi)
    nu_2 = sqrt(k_t / ((r_eq ** 2) * 1.0e-20 * mu_2)) / (2. * pi)

    print(nu_1 * 33.35e-12)
    print(nu_2 * 33.35e-12)

main()
