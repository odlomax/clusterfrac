from maths.fields.gaussian_random_field import scalar_grf
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import numpy as np
from maths.random.probability_density_function import pdf
#from maths.graphs.mst import euclidean_mst


grid_size=(500,500)
image_rescale=0.2

beta=(4.0,3.0,2.0)
sigma=(0.5,1.0,2.0)
seed=1

np.random.seed(seed)
A_11=scalar_grf(grid_size,3.)
A_11.normalise(1.,True)
shift=A_11.com_shift()



A_list=[]

for i in beta:
    for j in sigma:
        np.random.seed(seed)
        A=scalar_grf(grid_size,i)
        A.normalise(j,True)
        A.com_shift(shift)
        A_list.append(A)


fig_1=plt.figure(figsize=(8,8))
ax=[fig_1.add_subplot(3,3,i+1) for i in range(9)]


for i in range(len(ax)):
    
    ax[i].set_aspect("equal")
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    ax[i].imshow(A_list[i].signal.real.T**image_rescale,vmin=0.,origin="lower",cmap=plt.cm.inferno)

at=AnchoredText(r"$H=1.0, \sigma=0.5$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[0].add_artist(at)

at=AnchoredText(r"$H=1.0, \sigma=1.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[1].add_artist(at)

at=AnchoredText(r"$H=1.0, \sigma=2.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[2].add_artist(at)

at=AnchoredText(r"$H=0.5, \sigma=0.5$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[3].add_artist(at)

at=AnchoredText(r"$H=0.5, \sigma=1.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[4].add_artist(at)

at=AnchoredText(r"$H=0.5, \sigma=2.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[5].add_artist(at)

at=AnchoredText(r"$H=0.0, \sigma=0.5$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[6].add_artist(at)

at=AnchoredText(r"$H=0.0, \sigma=1.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[7].add_artist(at)

at=AnchoredText(r"$H=0.0, \sigma=2.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[8].add_artist(at)

fig_1.tight_layout()
fig_1.subplots_adjust(hspace=0,wspace=0)
fig_1.savefig("fbm.png",dpi=150)



fig_2=plt.figure(figsize=(8,8))
ax=[fig_2.add_subplot(3,3,i+1) for i in range(9)]


for i in range(len(ax)):
    
    A_pdf=pdf(A_list[i].signal.real)
    x=A_pdf.random(300)
    """
    mst=euclidean_mst(x)
    start,finish=mst.mst_edge_positions()
    """
    ax[i].set_xlim(0.,1.)
    ax[i].set_ylim(0.,1.)
    ax[i].set_aspect("equal")
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    """
    for j in range(start.shape[0]):
        ax[i].plot([start[j,0],finish[j,0]],[start[j,1],finish[j,1]],color="darkorange")
    """
    ax[i].plot(x[:,0],x[:,1],"o",markersize=2.,color="black")


    

at=AnchoredText(r"$H=1.0, \sigma=0.5$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[0].add_artist(at)

at=AnchoredText(r"$H=1.0, \sigma=1.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[1].add_artist(at)

at=AnchoredText(r"$H=1.0, \sigma=2.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[2].add_artist(at)

at=AnchoredText(r"$H=0.5, \sigma=0.5$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[3].add_artist(at)

at=AnchoredText(r"$H=0.5, \sigma=1.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[4].add_artist(at)

at=AnchoredText(r"$H=0.5, \sigma=2.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[5].add_artist(at)

at=AnchoredText(r"$H=0.0, \sigma=0.5$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[6].add_artist(at)

at=AnchoredText(r"$H=0.0, \sigma=1.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[7].add_artist(at)

at=AnchoredText(r"$H=0.0, \sigma=2.0$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[8].add_artist(at)


fig_2.tight_layout()
fig_2.subplots_adjust(hspace=0,wspace=0)
fig_2.savefig("cluster.png",dpi=150)












