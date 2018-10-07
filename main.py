from sys import argv
from os import listdir
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def get_theta (filename) :
  theta = filename.split("theta")[1][:-4]
  return float(theta)

def get_r (filename) :
  r = filename.split("theta")[0].split("r")[-1]
  return float(r)

def get_energy (filename) :
  f = open(dir + "/" +filename, "r")
  for line in f :
    if "SCF Done" in line :
      l = line.split()
      return float(l[4])

dir = argv[1]
file_list = listdir(dir)

r      = [get_r(filename)      for filename in file_list]
theta  = [get_theta(filename)  for filename in file_list]
energy = [get_energy(filename) for filename in file_list] 

ax = plt.axes(projection = "3d")
ax.plot_trisurf( r, theta, energy, cmap = "inferno" )
plt.show()
