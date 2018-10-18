import sys
import getopt
from os import listdir
from math import *
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np


def main():
    # Read in command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:r:t:")
    except getopt.GetoptError as err:
        print(err + "\n")
        help()
        sys.exit()

    # Check that a directory was supplied
    if len(args) == 0:
        print("No directory supplied.\n")
        help()
        sys.exit()

    # Read in data
    dir = args[0]
    data = get_data(dir)

    # Default values for options
    r_tol = 0.11
    t_tol = 20.
    out   = None

    # Process supplied arguments for options
    for o, a in opts:
        if o == '-h':
            help()
            sys.exit()
        elif o == '-r':
            try:
                r_tol = float(a)
            except ValueError:
                print("r_tol must be a number. Using default.")
        elif o == '-t':
            try:
                t_tol = float(a)
            except ValueError:
                print("t_tol must be a number. Using default.")
        elif o == '-n':
            out = a

    # Generate filename if none supplied
    if out == None:
        print("No filename supplied. Using default.")
        out = "%s_r%.2f_t%.1f" % (dir.split('/')[-1], r_tol, t_tol)
    print("Proceeding with r_tol = %f and t_tol = %f.\n" % (r_tol, t_tol))

    # Get equilibrium geometry
    r_eq, t_eq, e_eq = min(data, key = lambda tup: tup[2])

    # Constants used in calculations
    wavenumber = 3.335641e-11  # Convert Hz to cm-1
    angstrom = 1.0e-10         # Convert angstrom to metres
    degree = pi / 180.         # Convert degrees to radians
    hartree = 4.35974e-18      # Convert hartree to J
    mu_1 = 2.0 * 1.66e-27      # Reduced mass for streching
    mu_2 = 0.5 * 1.66e-27      # Reduced mass for bending

    # Fit data along symmetric stretch
    r_mode = [(r, e) for (r, t, e) in data if t == t_eq]  # Keep t = t_eq)
    k_r = 2. * fit_mode(r_mode, r_tol, r_eq, "STRETCH_%s.png" % out,
        "r / Angstrom", "Energy / Hartrees")
    # k in hartrees angstrom^-2. Convert to J m^-2
    k_r *= hartree / ((angstrom)**2)

    # Fit data along bending mode
    t_mode = [(t, e) for (r, t, e) in data if r == r_eq]  # Keep r = r_eq
    k_t = 2. * fit_mode(t_mode, t_tol, t_eq, "BEND_%s.png" % out,
        "theta / Degrees", "Energy / Hartrees")
    # k in hartrees degree^-2. Convert to J radian^-2
    k_t *= hartree / (degree**2)

    # Calculate vibrational frequencties in wavenumbers
    nu_1 = sqrt(k_r / mu_1) / (2. * pi)  * wavenumber
    nu_2 = sqrt(k_t / ((r_eq * angstrom)**2 * mu_2)) / (2. * pi) * wavenumber

    write_output(r_eq, t_eq, e_eq, nu_1, nu_2, r_tol, t_tol, dir,
        "OUTPUT_%s.txt" % out)
    plot_surface(data, "SURFACE_%s.png" % out)


def get_data(dir):
    try:
        file_list = [str for str in listdir(dir) if str.endswith('.out')]
    except (FileNotFoundError, NotADirectoryError):
        print("Directory %s not found." % dir)
        sys.exit()

    data = [get_point(filename, dir) for filename
        in file_list if get_point(filename, dir) != None]

    if len(data) != 0:
        print("Loaded %d files \n" % len(data))
        return data
    else:
        print("No .out files found in given directory")
        sys.exit()


def get_point(filename, dir):
    try:
        r = get_r(filename)
        t = get_theta(filename)
        e = get_energy(dir + '/' + filename)
        return (r,t,e)
    except:
        print("Failed to get data point for " + filename)


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


def write_output(r_eq, t_eq, e_eq, nu_1, nu_2, r_tol, t_tol, dir, out):
    # Output data to file
    f = open(out, 'w')

    f.write("Data read from %s \n\n\n" % dir)
    f.write("Equilibrium geometry: \n\n")
    f.write("Bond length = %.2f Angstroms\n" % r_eq)
    f.write("Bond angle = %.0f degrees\n" % t_eq)
    f.write("Energy = %f Hartrees\n\n\n" % e_eq)
    f.write("Vibrational frequencties:\n\n")
    f.write("Frequency of symmetric stretch = %f cm^-1\n" % nu_1)
    f.write("An r tolerance of %.2f Angstrom was used.\n\n" % r_tol)
    f.write("Frequency of bending mode = %f cm^-1\n" % nu_2)
    f.write("A theta tolerance of %.1f degrees was used." % t_tol)
    f.close()
    print("Wrote calculated parameters data to " + out)


def plot_surface(data, out):
    # Prepare points to be plotted
    # Values are interpolated to give smoother colour gradient on plot
    r,t,e = zip(*data)
    r_space = np.linspace(min(r),max(r),100)
    t_space = np.linspace(min(t),max(t),100)
    grid_r, grid_t = np.meshgrid(r_space, t_space)
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
    plt.savefig(out, dpi = 300)
    print("Saved surface plot as " + out)
    plt.close()


def fit_mode(data, tol, eq, out, xlabel, ylabel):
    data.sort(key = lambda tup: tup[0])
    ys = [y for (x,y) in data if abs(x - eq) < tol]
    xs = [x for (x,y) in data if abs(x - eq) < tol]
    p = np.polyfit(xs, ys, 2)

    xs_fit = np.linspace(min(xs), max(xs), 100)
    ys_fit = [np.polyval(p, x) for x in xs_fit]

    plt.scatter(xs, ys, marker='+', linewidth=0.5,color='k')
    plt.plot(xs_fit, ys_fit, linewidth=0.8, color='r')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(out, dpi = 300)

    print("Saved plot of fitted mode to " + out)
    plt.close()
    return p[0]


def help():
    print("Directory containing .out files must be supplied as argument after "
        "all options")
    print("-r NUMBER : (OPTIONAL) set tolerance in r (in Angstroms) value "
        "when making quadratic approximation for fitting. Default = 0.11")
    print("-t NUMBER : (OPTIONAL) set tolerance in theta (in degrees) value "
        "when making quadratic approximation for fitting. Default = 20.0")
    print("-n STRING : (OPTIONAL) set filename for output files. Default "
        "from other parameters if not supplied")
    print("-h : List command line options")


if __name__ == "__main__":
    main()
