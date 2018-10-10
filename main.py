from sys import argv
from os import listdir

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
        print("Failed to get data point for " + filename)


#Functions for processing data
def fit_data(data):
    r_eq, t_eq, e_eq = min(data, key = lambda tup: tup[2])

    r_fit = [(r, t, e) for (r, t, e) in data if r == r_eq]
    r_fit.sort(key = lambda tup: tup[0])

    t_fit = [(r, t, e) for (r, t, e) in data if t == t_eq]
    t_fit.sort(key = lambda tup: tup[1])

main()
