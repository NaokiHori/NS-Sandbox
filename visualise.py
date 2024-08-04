import os
import numpy as np
from matplotlib import pyplot

root = "output/save"

dnames = [f"{root}/{dname}" for dname in os.listdir(root)]
dnames = [dname for dname in dnames if os.path.isdir(dname)]
dnames = sorted(dnames)

fig = pyplot.figure(figsize=(8, 8))
axes = [fig.add_subplot(211), fig.add_subplot(212)]
for dname in dnames:
    ux = np.load(f"{dname}/ux.npy")
    uy = np.load(f"{dname}/uy.npy")
    axes[0].clear()
    axes[1].clear()
    axes[0].contourf(ux)
    axes[1].contourf(uy)
    keywords = {
            "aspect": "equal",
            }
    axes[0].set(**keywords)
    axes[1].set(**keywords)
    pyplot.show(block=False)
    pyplot.pause(1.e-1)
pyplot.close()

