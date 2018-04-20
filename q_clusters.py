import numpy as np
import pandas as pd
from maths.points.box_fractal import gw04_box_fractal
from maths.points.radial_profile import power_law_profile
from clusterfrac.cluster import star_cluster
from maths.tensors.rotation_matrix import r_matrix_3d
from matplotlib import pyplot as plt


n_dim=3
n_run=10

n_star_min=316
n_star_max=1000

alpha_min=0.
alpha_max=3.

D_min=1.
D_max=3.




n_star=np.random.uniform(n_star_min,n_star_max,n_run).astype(np.int)
alpha=np.random.uniform(alpha_min,alpha_max,n_run)
D=np.random.uniform(D_min,D_max,n_run)

# initialise observables
R_hull=np.zeros(n_run)
m_mean=np.zeros(n_run)
m_m2=np.zeros(n_run)
m_m3=np.zeros(n_run)
m_m4=np.zeros(n_run)
s_mean=np.zeros(n_run)
s_m2=np.zeros(n_run)
s_m3=np.zeros(n_run)
s_m4=np.zeros(n_run)


for i in range(n_run):
    
    print("fractal cluster",i+1,"of",n_run)
    
    N_gen=np.ceil(np.log2(n_star[i])/D[i]).astype(np.int)
    
    x=gw04_box_fractal(n_dim,D[i],N_gen,True)
    
    indices=np.arange(x.shape[0])
    in_sphere=(x**2).sum(axis=1)<1.
    #x=x[indices[in_sphere],:]
    
    n_star[i]=x.shape[0]
    
    
    phi=np.random.uniform(0.,2.*np.pi)
    theta=np.arccos(np.random.uniform(-1.,1.))
    
    R_1=r_matrix_3d([0.,0.,1.],-phi)
    R_2=r_matrix_3d([0.,1.,0.],-theta)

    x=(R_1.rotate_vec(x.T)).T
    x=(R_2.rotate_vec(x.T)).T
    
    
    cluster=star_cluster(x[:,:2])
    R_hull[i]=cluster.hull_length()
    m_mean[i]=cluster.m_bar()
    m_m2[i]=cluster.m_moment(2)
    m_m3[i]=cluster.m_moment(3)
    m_m4[i]=cluster.m_moment(4)
    s_mean[i]=cluster.s_bar()
    s_m2[i]=cluster.s_moment(2)
    s_m3[i]=cluster.s_moment(3)
    s_m4[i]=cluster.s_moment(4)
    
    plt.figure(i)
    plt.axes().set_xlim(-1,1)
    plt.axes().set_ylim(-1,1)
    plt.axes().set_aspect("equal")
    plt.axes().set_xticks([])
    plt.axes().set_yticks([])
    plt.plot(x[:,0],x[:,1],".")
    plt.tight_layout()
    plt.show()
    
    
# make table of data
data_array=np.array((n_star,D,R_hull,m_mean,m_m2,m_m3,m_m4,s_mean,s_m2,s_m3,s_m4)).T
header=("N_star","D","R_hull","m_mean","m_m2","m_m3","m_m4","s_mean","s_m2","s_m3","s_m4")

fractal_data=pd.DataFrame(data=data_array,columns=header)
#fractal_data.to_csv("fractal_cluster.dat") 


for i in range(n_run):
    
    print("radial cluster",i+1,"of",n_run)
    x=power_law_profile(n_dim,alpha[i],n_star[i])
    cluster=star_cluster(x[:,:2])
    R_hull[i]=cluster.hull_length()
    m_mean[i]=cluster.m_bar()
    m_m2[i]=cluster.m_moment(2)
    m_m3[i]=cluster.m_moment(3)
    m_m4[i]=cluster.m_moment(4)
    s_mean[i]=cluster.s_bar()
    s_m2[i]=cluster.s_moment(2)
    s_m3[i]=cluster.s_moment(3)
    s_m4[i]=cluster.s_moment(4)
    
    plt.figure(i+n_run)
    plt.axes().set_xlim(-1,1)
    plt.axes().set_ylim(-1,1)
    plt.axes().set_aspect("equal")
    plt.axes().set_xticks([])
    plt.axes().set_yticks([])
    plt.plot(x[:,0],x[:,1],".")
    plt.tight_layout()
    plt.show()
    
# make table of data
data_array=np.array((n_star,alpha,R_hull,m_mean,m_m2,m_m3,m_m4,s_mean,s_m2,s_m3,s_m4)).T
header=("N_star","alpha","R_hull","m_mean","m_m2","m_m3","m_m4","s_mean","s_m2","s_m3","s_m4")

radial_data=pd.DataFrame(data=data_array,columns=header)
#radial_data.to_csv("radial_cluster.dat")
