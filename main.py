from sys import argv
from os import listdir
from datetime import datetime
from math import *
from numpy import *
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm


def main():
    try:
        dir = argv[1]
    except IndexError:
        print("Please enter directory as command line "
            "containing output files as argument.")
        return # Exit program

    file_list = [str for str in listdir(dir) if str.endswith('.out')]

    if len(file_list) == 0:
        print("No .out files found in given directory")
        return # Exit program

    data = [get_data(filename, dir) for filename
        in file_list if get_data(filename, dir) != None]

    print("Loaded " + str(len(data)) + " files")

    dt = datetime.now()
    out = "%d-%02d-%02d_%02d:%02d:%02d" % (dt.year,
                                           dt.month,
                                           dt.day,
                                           dt.hour,
                                           dt.minute,
                                           dt.second)
    plot(data, "FIGURE_" + out + ".png")

    params = fit_data(data, .95, .99)

    f = open("OUTPUT_" + out + ".txt", 'w')

    f.write("Data read from %s \n\n\n" % dir)
    f.write("Equilibrium geometry: \n\n")
    f.write("Bond length = %.2f Angstroms\n" % params.get('r_eq'))
    f.write("Bond angle = %.0f degrees\n" % params.get('t_eq'))
    f.write("Energy = %f Hartrees\n\n\n" % params.get('e_eq'))
    f.write("Vibrational frequencties:\n\n")
    f.write("Frequency of symmetric stretch = %f cm^-1\n"
        % params.get('nu_1'))
    f.write("Symmetric stretch data fitted using %d data points with R^2 = %f\n"
        % (params.get('r_fit')[2], params.get('r_fit')[1]))
    f.write("\n")
    f.write("Frequency of bending mode = %f cm^-1\n"
        % params.get('nu_2'))
    f.write(("Bending mode data fitted using %d data points with R^2 = %f\n")
        % (params.get('t_fit')[2], params.get('t_fit')[1]))
    f.close()
    print("Wrote calculated parameters data to OUTPUT_%s.txt" % out)



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
        if "SCF Done" in line:
            l = line.split()
            f.close()
            return float(l[4])


def get_data(filename, dir):
    try:
        r = get_r(filename)
        t = get_theta(filename)
        e = get_energy(dir + '/' + filename)
        return (r,t,e)
    except:
        print("Failed to get data point for " + filename)


# Plots data on surface plots and outputs png file of figure
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
    ax.set_xlabel("r / Angstroms")
    ax.set_ylabel("theta / degrees")
    ax.set_zlabel("Energy / Hartrees")
    ax.grid(False)
    ax.plot_surface(grid_r,grid_t,grid_e, rstride=1, cstride=1,
        cmap='inferno')
    plt.savefig(out, dpi = 200)
    print("Saved surface plot as " + out)


# Functions for fitting data

def fit_data(data, tol_r, tol_t):
    # Get equilibrium geometry
    r_eq, t_eq, e_eq = min(data, key = lambda tup: tup[2])

    # Fit data along symmetric stretch
    r_mode = [(r, t, e) for (r, t, e) in data if t == t_eq ]  # Keep t = t_eq
    r_mode.sort(key = lambda tup: tup[0])
    rs = [x[0] * 1.0e-10 for x in r_mode]      # Get r in m
    es = [x[2] * 4.35974e-18 for x in r_mode]  # Get E in J
    r_fit = fit_to_tol(rs, es, r_eq * 1.0e-10, tol_r)
    k_r = 2. * r_fit[0]

    # Fit data along bending mode
    t_mode = [(r, t, e) for (r, t, e) in data if r == r_eq]  # Keep r = r_eq
    t_mode.sort(key = lambda tup: tup[1])
    ts = [x[1] * pi / 180. for x in t_mode]    # Get theta in radians
    es = [x[2] * 4.35974e-18 for x in t_mode]  # Get E in J
    t_fit = fit_to_tol(ts, es, t_eq * pi / 180., tol_t)
    k_t = 2. * t_fit[0]

    # Calculate vibrational frequencties in wavenumbers
    mu_1 = 2.0 * 1.66e-27
    mu_2 = 0.5 * 1.66e-27

    nu_1 = sqrt(k_r / mu_1) / (2. * pi) * 33.35641e-12
    nu_2 = sqrt(k_t / ((r_eq ** 2) * 1.0e-20 * mu_2)) / (2. * pi) * 33.35641e-12

    return {'r_eq' : r_eq,
            't_eq' : t_eq,
            'e_eq' : e_eq,
            'r_fit': r_fit,
            't_fit': t_fit,
            'k_r'  : k_r,
            'k_t'  : k_t,
            'nu_1' : nu_1,
            'nu_2' : nu_2
            }


# Fits quadratic to xs and ys, reducing maximum displacement from equilibrium
# Until R^2 rises above minimum value
def fit_to_tol(xs, ys, x_eq, tol):
    p = polyfit(xs, ys, 2) # Fit data to quadratic
    range = max([abs(x-x_eq) for x in xs])
    step = range / (2. * len(xs)) # Amount max displacement reduced
                                  # by after each fit
    R2 = get_r2(p, xs, ys)

    while R2 < tol:     # R^2 below threshold
        #Reduce maximum displacement from equilibrium
        range -= step
        ys = [y for (i, y) in enumerate(ys) if abs(xs[i] - x_eq) < range]
        xs = [x for x in xs if abs(x - x_eq) < range]

        #Recalculate fit
        p = polyfit(xs, ys, 2)
        R2 = get_r2(p,xs,ys)

    return (p[0], R2, len(xs))


# Returns r^2 value for polynomial p which is fitted to xs and ys
def get_r2(p, xs, ys):
    y_fit = [polyval(p, x) for x in xs] # y values predicted by parameters
    y_avg = sum(ys)/len(ys)             # Average y value
    SStot = sum((ys - y_avg)**2)        # Sum of squares
    y_res = subtract(ys, y_fit)         # Residuals
    SSres = sum(y_res**2)               # Residual sum of squares

    return 1. - SSres / SStot

main()
