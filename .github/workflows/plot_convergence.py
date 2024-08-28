import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

root = sys.argv[1]
schemes = ["advx", "advy", "difx", "dify", "pres"]

fig = pyplot.figure()
ax = fig.add_subplot()

for cnt, scheme in enumerate(schemes):
    file_name = f"{root}/{scheme}.dat"
    data = np.loadtxt(file_name)
    if 0 == cnt:
        ax.plot(data[:, 0], data[:, 0]**-2., label="2nd order")
    ax.plot(data[:, 0], data[:, 1], label=f"{scheme}, L^2",   marker="o")
    ax.plot(data[:, 0], data[:, 2], label=f"{scheme}, L^inf", marker="x")

ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

pyplot.savefig(sys.argv[2], dpi=300)
pyplot.close()

