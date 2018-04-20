import numpy as np
from maths.graphs.mst import euclidean_mst
from matplotlib import pyplot as plt


x=np.random.random((30,2))


mst=euclidean_mst(x)

s_0,s_1=mst.com_edge_positions()
m_0,m_1=mst.mst_edge_positions()



plt.figure(1,figsize=(6,6))

plt.axes().set_xlim(0.,1.)
plt.axes().set_ylim(0.,1.)
plt.axes().set_aspect("equal")
plt.axes().set_xticks([])
plt.axes().set_yticks([])

plt.plot(x[:,0],x[:,1],".")

plt.tight_layout()
plt.show()

plt.figure(2,figsize=(6,6))

plt.axes().set_xlim(0.,1.)
plt.axes().set_ylim(0.,1.)
plt.axes().set_aspect("equal")
plt.axes().set_xticks([])
plt.axes().set_yticks([])

for i in range(m_0.shape[0]):
    plt.plot([m_0[i,0],m_1[i,0]],[m_0[i,1],m_1[i,1]],color="black")

plt.plot(x[:,0],x[:,1],".")

plt.tight_layout()
plt.show()

plt.figure(3,figsize=(6,6))

plt.axes().set_xlim(0.,1.)
plt.axes().set_ylim(0.,1.)
plt.axes().set_aspect("equal")
plt.axes().set_xticks([])
plt.axes().set_yticks([])

for i in range(s_0.shape[0]):
    plt.plot([s_0[i,0],s_1[i,0]],[s_0[i,1],s_1[i,1]],color="black")

plt.plot(x[:,0],x[:,1],".")

plt.tight_layout()
plt.show()