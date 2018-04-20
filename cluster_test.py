from clusterfrac.cluster import star_cluster
from clusterfrac.model import cluster_model
from clusterfrac.estimator import param_estimator
from maths.points.ra_dec import ra_dec_project
from maths.points.fuse_points import fuse_close_companions
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse

"""
model_2d2d=cluster_model("2d2d")


data_2d2d_0=model_2d2d.make_training_data(32,100,15000)
data_2d2d_0.to_csv("data_2d2d_0.dat")


data_2d2d_1=model_2d2d.make_training_data(100,316,15000)
data_2d2d_1.to_csv("data_2d2d_1.dat")


data_2d2d_2=model_2d2d.make_training_data(316,1000,15000)
data_2d2d_2.to_csv("data_2d2d_2.dat")


model_3d2d=cluster_model("3d2d")
data_3d2d_0=model_3d2d.make_training_data(32,100,15000)
data_3d2d_0.to_csv("data_3d2d_0.dat")

data_3d2d_1=model_3d2d.make_training_data(100,316,15000)
data_3d2d_1.to_csv("data_3d2d_1.dat")

data_3d2d_2=model_3d2d.make_training_data(316,1000,15000)
data_3d2d_2.to_csv("data_3d2d_2.dat")


model_3d3d=cluster_model("3d3d")
data_3d3d_0=model_3d3d.make_training_data(32,100,15000)
data_3d3d_0.to_csv("data_3d3d_0.dat")

data_3d3d_1=model_3d3d.make_training_data(100,316,15000)
data_3d3d_1.to_csv("data_3d3d_1.dat")

data_3d3d_2=model_3d3d.make_training_data(316,1000,15000)
data_3d3d_2.to_csv("data_3d3d_2.dat")
"""

data_2d2d_0=pd.read_csv("data_2d2d_0.dat")
estimator_2d2d_0=param_estimator(data_2d2d_0)

data_2d2d_1=pd.read_csv("data_2d2d_1.dat")
estimator_2d2d_1=param_estimator(data_2d2d_1)

data_2d2d_2=pd.read_csv("data_2d2d_2.dat")
estimator_2d2d_2=param_estimator(data_2d2d_2)


data_3d2d_0=pd.read_csv("data_3d2d_0.dat")
estimator_3d2d_0=param_estimator(data_3d2d_0)

data_3d2d_1=pd.read_csv("data_3d2d_1.dat")
estimator_3d2d_1=param_estimator(data_3d2d_1)

data_3d2d_2=pd.read_csv("data_3d2d_2.dat")
estimator_3d2d_2=param_estimator(data_3d2d_2)

radial_cluster=pd.read_csv("radial_cluster.dat")
fractal_cluster=pd.read_csv("fractal_cluster.dat")


taurus_dist=140
chai_dist=160
ic348_dist=315
lupus3_dist=170
ophiuchus_dist=130
ic2391_dist=150

cc_dist=0.005





taurus_pos=np.loadtxt("yso_taurus.txt")
taurus_x=ra_dec_project(taurus_pos[:,:3],taurus_pos[:,3:])[:,:2]*taurus_dist
taurus_x=fuse_close_companions(taurus_x,cc_dist)
taurus=star_cluster(taurus_x)

chai_pos=np.loadtxt("yso_chai.txt")
chai_x=ra_dec_project(chai_pos[:,:3],chai_pos[:,3:])[:,:2]*chai_dist
chai_x=fuse_close_companions(chai_x,cc_dist)
chai=star_cluster(chai_x)

ic348_pos=np.loadtxt("yso_ic348.txt")
ic348_x=ra_dec_project(ic348_pos[:,:3],ic348_pos[:,3:])[:,:2]*ic348_dist
ic348_x=fuse_close_companions(ic348_x,cc_dist)
ic348=star_cluster(ic348_x)

lupus3_pos=np.loadtxt("yso_lupus3.txt")
lupus3_x=ra_dec_project(lupus3_pos[:,:3],lupus3_pos[:,3:])[:,:2]*lupus3_dist
lupus3_x=fuse_close_companions(lupus3_x,cc_dist)
lupus3=star_cluster(lupus3_x)

ophiuchus_pos=np.loadtxt("yso_ophiuchus.txt")
ophiuchus_x=ra_dec_project(ophiuchus_pos[:,:3],ophiuchus_pos[:,3:])[:,:2]*ophiuchus_dist
ophiuchus_x=fuse_close_companions(ophiuchus_x,cc_dist)
ophiuchus=star_cluster(ophiuchus_x)

ic2391_pos=np.loadtxt("stars_ic2391.txt")
ic2391_x=ra_dec_project(ic2391_pos[:,:3],ic2391_pos[:,3:])[:,:2]*ic2391_dist


ic2391_x=fuse_close_companions(ic2391_x,cc_dist)


ic2391=star_cluster(ic2391_x)

taurus_table=taurus.make_table_row()
chai_table=chai.make_table_row()
ic348_table=ic348.make_table_row()
lupus3_table=lupus3.make_table_row()
ophiuchus_table=ophiuchus.make_table_row()
ic2391_table=ic2391.make_table_row()


estimator_3d2d_1.estimate_params(taurus_table)

print("taurus")
print(taurus_table)
print()


estimator_3d2d_1.estimate_params(chai_table)

print("chai")
print(chai_table)
print()

estimator_3d2d_2.estimate_params(ic348_table)

print("ic348")
print(ic348_table)
print()

estimator_3d2d_0.estimate_params(lupus3_table)

print("lupus3")
print(lupus3_table)

observations=np.zeros((6,2))

estimator_3d2d_1.estimate_params(ophiuchus_table)

print("ophiuchus")
print(ophiuchus_table)
print()

estimator_3d2d_1.estimate_params(ic2391_table)

print("ic2391")
print(ic2391_table)
print()

observations[5,:]=[taurus_table["s_mean"]/taurus_table["R_hull"],taurus_table["m_mean"]/taurus_table["R_hull"]]
observations[4,:]=[chai_table["s_mean"]/chai_table["R_hull"],chai_table["m_mean"]/chai_table["R_hull"]]
observations[1,:]=[ic348_table["s_mean"]/ic348_table["R_hull"],ic348_table["m_mean"]/ic348_table["R_hull"]]
observations[0,:]=[lupus3_table["s_mean"]/lupus3_table["R_hull"],lupus3_table["m_mean"]/lupus3_table["R_hull"]]
observations[2,:]=[ophiuchus_table["s_mean"]/ophiuchus_table["R_hull"],ophiuchus_table["m_mean"]/ophiuchus_table["R_hull"]]
observations[3,:]=[ic2391_table["s_mean"]/ic2391_table["R_hull"],ic2391_table["m_mean"]/ic2391_table["R_hull"]]



"""
fig_1,ax_1=plt.subplots()

D_ax=ax_1.scatter(np.sqrt(np.pi)*fractal_cluster["s_mean"]/fractal_cluster["R_hull"],fractal_cluster["m_mean"]/fractal_cluster["R_hull"],marker=".",c=fractal_cluster["D"],cmap=plt.cm.viridis)
alpha_ax=ax_1.scatter(np.sqrt(np.pi)*radial_cluster["s_mean"]/radial_cluster["R_hull"],radial_cluster["m_mean"]/radial_cluster["R_hull"],marker=".",c=radial_cluster["alpha"],cmap=plt.cm.inferno)
bbox_props=dict(boxstyle="circle",fc="white",pad=0.2)

for i in range(observations.shape[0]):
    ax_1.text(np.sqrt(np.pi)*observations[i,0],observations[i,1],i+1,ha="center",va="center",bbox=bbox_props)

ax_1.set_xlim(0.15,1.25)
ax_1.set_ylim(0.15,0.75)
ax_1.set_xlabel(r"$\overline{s}$")
ax_1.set_ylabel(r"$\overline{m}$")



D_ax.set_clim(1.,3.)
alpha_ax.set_clim(0.,3.)

D_cb=fig_1.colorbar(D_ax)
D_cb.set_label(r"$D$")
D_cb.set_ticks(np.linspace(1.,3.,6))


alpha_cb=fig_1.colorbar(alpha_ax)
alpha_cb.set_label(r"$\alpha$")
alpha_cb.set_ticks(np.linspace(0,3,6))


fig_1.tight_layout()



fig_1.savefig("alpha_D_plot.png",dpi=150)
"""



fig_2=plt.figure(figsize=(6.,12.))
ax=[fig_2.add_subplot(3,1,i+1) for i in range(3)]


D_ax=ax[0].scatter(np.sqrt(np.pi)*fractal_cluster["s_mean"]/fractal_cluster["R_hull"],fractal_cluster["m_mean"]/fractal_cluster["R_hull"],marker=".",c=fractal_cluster["D"],cmap=plt.cm.viridis)
alpha_ax=ax[0].scatter(np.sqrt(np.pi)*radial_cluster["s_mean"]/radial_cluster["R_hull"],radial_cluster["m_mean"]/radial_cluster["R_hull"],marker=".",c=radial_cluster["alpha"],cmap=plt.cm.inferno)

bbox_props=dict(boxstyle="circle",fc="white",pad=0.2)

for i in range(observations.shape[0]):
    ax[0].text(np.sqrt(np.pi)*observations[i,0],observations[i,1],i+1,ha="center",va="center",bbox=bbox_props)
    
ax[0].set_xlim(0.15,1.25)
ax[0].set_ylim(0.15,0.75)
ax[0].set_ylabel(r"$\overline{m}$")
ax[0].set_xticks([])

D_ax.set_clim(1.,3.)
alpha_ax.set_clim(0.,3.)


H_ax=ax[1].scatter(np.sqrt(np.pi)*data_2d2d_1["s_mean"]/data_2d2d_1["R_hull"],data_2d2d_1["m_mean"]/data_2d2d_1["R_hull"],marker=".",c=data_2d2d_1["H"],cmap=plt.cm.inferno)

for i in range(observations.shape[0]):
    ax[1].text(np.sqrt(np.pi)*observations[i,0],observations[i,1],i+1,ha="center",va="center",bbox=bbox_props)
    
ax[1].set_xlim(0.15,1.25)
ax[1].set_ylim(0.15,0.75)
ax[1].set_ylabel(r"$\overline{m}$")
ax[1].set_xticks([])

H_ax.set_clim(0.,1.)

sigma_ax=ax[2].scatter(np.sqrt(np.pi)*data_2d2d_1["s_mean"]/data_2d2d_1["R_hull"],data_2d2d_1["m_mean"]/data_2d2d_1["R_hull"],marker=".",c=data_2d2d_1["sigma"],cmap=plt.cm.inferno)

for i in range(observations.shape[0]):
    ax[2].text(np.sqrt(np.pi)*observations[i,0],observations[i,1],i+1,ha="center",va="center",bbox=bbox_props)


ax[2].set_xlim(0.15,1.25)
ax[2].set_ylim(0.15,0.75)
ax[2].set_xlabel(r"$\overline{s}$")
ax[2].set_ylabel(r"$\overline{m}$")

sigma_ax.set_clim(0.5,2.5)


fig_2.tight_layout()
fig_2.subplots_adjust(hspace=0,right=0.85)


pos=ax[0].get_position()
cbar_ax = fig_2.add_axes([0.87, pos.y0, 0.04, 0.5*(pos.y1-pos.y0)])
D_cb=fig_2.colorbar(D_ax, cax=cbar_ax,ticks=np.linspace(1.2,2.8,5),label=r"$D$")

pos=ax[0].get_position()
cbar_ax = fig_2.add_axes([0.87, pos.y0+0.5*(pos.y1-pos.y0), 0.04, 0.5*(pos.y1-pos.y0)])
alpha_cb=fig_2.colorbar(alpha_ax, cax=cbar_ax,ticks=np.linspace(0.3,2.7,5),label=r"$\alpha$")


pos=ax[1].get_position()
cbar_ax = fig_2.add_axes([0.87, pos.y0, 0.04, pos.y1-pos.y0])
H_cb=fig_2.colorbar(H_ax, cax=cbar_ax,ticks=np.linspace(0.1,0.9,5),label=r"$H$")


pos=ax[2].get_position()
cbar_ax = fig_2.add_axes([0.87, pos.y0, 0.04, pos.y1-pos.y0])
sigma_cb=fig_2.colorbar(sigma_ax, cax=cbar_ax,ticks=np.linspace(0.7,2.3,5),label=r"$\sigma$")




fig_2.savefig("H_sigma_plot.png",dpi=150)

"""

fig_3=plt.figure(figsize=(12.,6.))
ax=[fig_3.add_subplot(1,3,i+1) for i in range(3)]


"""

fig_3,ax=plt.subplots(1,4,figsize=(12.,6.),gridspec_kw = {'width_ratios':[1, 1,2,2]})

ax[1].set_yticks([])
ax[2].set_yticks([])
ax[3].set_yticks([])

D_ax=ax[0].scatter(fractal_cluster["D"],fractal_cluster["m_mean"]/(np.sqrt(np.pi)*fractal_cluster["s_mean"]),c="indigo",marker=".")

ax[0].set_xlim(1,3)
ax[0].set_ylim(0,2)
ax[0].set_xticks(np.linspace(1.2,2.8,5))
ax[0].set_xlabel(r"$D$",)
ax[0].set_ylabel(r"$Q$")
ax[0].set_yticks(np.linspace(0,2,6))
ax[0].spines["right"].set_visible(False)

for i in range(observations.shape[0]):
    ax[0].plot(np.linspace(1,3,10),[observations[i,1]/(np.sqrt(np.pi)*observations[i,0])]*10,color="teal",linewidth=2)



alpha_ax=ax[1].scatter((radial_cluster["alpha"]),radial_cluster["m_mean"]/(np.sqrt(np.pi)*radial_cluster["s_mean"]),c="orange",marker=".")

ax[1].set_xlim(0,3)
ax[1].set_ylim(0,2)
ax[1].set_xticks([])
ax[1].spines["left"].set_visible(False)

for i in range(observations.shape[0]):
    ax[1].plot(np.linspace(0,3,10),[observations[i,1]/(np.sqrt(np.pi)*observations[i,0])]*10,color="teal",linewidth=2)


ax_top=ax[1].twiny()
ax_top.set_xlim(0,3)
ax_top.set_xticks(np.linspace(0.3,2.7,5))
ax_top.set_xlabel(r"$\alpha$")
ax_top.spines["left"].set_visible(False)

H_ax=ax[2].scatter((data_2d2d_1["H"]),data_2d2d_1["m_mean"]/(np.sqrt(np.pi)*data_2d2d_1["s_mean"]),c=data_2d2d_1["sigma"],cmap=plt.cm.inferno,marker=".")

ax[2].set_xlim(0,1)
ax[2].set_ylim(0,2)
ax[2].set_xticks(np.linspace(0.1,0.9,5))
ax[2].set_xlabel(r"$H$")
H_ax.set_clim(0.5,2.5)



sigma_ax=ax[3].scatter((data_2d2d_1["sigma"]),data_2d2d_1["m_mean"]/(np.sqrt(np.pi)*data_2d2d_1["s_mean"]),c=data_2d2d_1["H"],cmap=plt.cm.inferno,marker=".")

ax[3].set_xlim(0.5,2.5)
ax[3].set_ylim(0,2)
sigma_ax.set_clim(0,1)
ax[3].set_xticks(np.linspace(0.7,2.3,5))
ax[3].set_xlabel(r"$\sigma$")



fig_3.tight_layout()
fig_3.subplots_adjust(wspace=0,top=0.85)


pos=ax[2].get_position()
cbar_ax = fig_3.add_axes([pos.x0, 0.87, pos.x1-pos.x0, 0.05])
H_cb=fig_3.colorbar(H_ax, cax=cbar_ax,ticks=np.linspace(0.7,2.3,5),label=r"$\sigma$",orientation="horizontal")
cbar_ax.xaxis.set_ticks_position("top")
cbar_ax.xaxis.set_label_position("top")

pos=ax[3].get_position()
cbar_ax = fig_3.add_axes([pos.x0, 0.87, pos.x1-pos.x0, 0.05])
sigma_cb=fig_3.colorbar(sigma_ax, cax=cbar_ax,ticks=np.linspace(0.1,0.9,5),label=r"$H$",orientation="horizontal")
cbar_ax.xaxis.set_ticks_position("top")
cbar_ax.xaxis.set_label_position("top")


fig_3.savefig("Q_plot.png",dpi=150)



fig_4,ax=plt.subplots(figsize=(6.,6.))


w,v=np.linalg.eig(estimator_2d2d_2.covar)

ell_size=2.*np.sqrt(w)

ell_angle=np.arctan2(v[1,0],v[0,0])*180./np.pi




ax.set_xlim(0,1)
ax.set_ylim(0.5,2.5)



ell=Ellipse((taurus_table["H_est"],taurus_table["sigma_est"]),ell_size[0],ell_size[1],ell_angle)
ell.set_alpha(0.3)

ax.add_artist(ell)
ax.text(taurus_table["H_est"],taurus_table["sigma_est"],6,ha="center",va="center",bbox=bbox_props)


ell=Ellipse((ic348_table["H_est"],ic348_table["sigma_est"]),ell_size[0],ell_size[1],ell_angle)
ell.set_alpha(0.3)

ax.add_artist(ell)
ax.text(ic348_table["H_est"],ic348_table["sigma_est"],2,ha="center",va="center",bbox=bbox_props)


w,v=np.linalg.eig(estimator_2d2d_1.covar)

ell_size=2.*np.sqrt(w)

ell_angle=np.arctan2(v[1,0],v[0,0])*180./np.pi

print(ell_size)



ell=Ellipse((chai_table["H_est"],chai_table["sigma_est"]),ell_size[0],ell_size[1],ell_angle)
ell.set_alpha(0.3)

ax.add_artist(ell)
ax.text(chai_table["H_est"],chai_table["sigma_est"],5,ha="center",va="center",bbox=bbox_props)


ell=Ellipse((ophiuchus_table["H_est"],ophiuchus_table["sigma_est"]),ell_size[0],ell_size[1],ell_angle)
ell.set_alpha(0.3)

ax.add_artist(ell)
ax.text(ophiuchus_table["H_est"],ophiuchus_table["sigma_est"],3,ha="center",va="center",bbox=bbox_props)

ell=Ellipse((ic2391_table["H_est"],ic2391_table["sigma_est"]),ell_size[0],ell_size[1],ell_angle)
ell.set_alpha(0.3)

ax.add_artist(ell)
ax.text(ic2391_table["H_est"],ic2391_table["sigma_est"],4,ha="center",va="center",bbox=bbox_props)


w,v=np.linalg.eig(estimator_2d2d_0.covar)

ell_size=2.*np.sqrt(w)

ell_angle=np.arctan2(v[1,0],v[0,0])*180./np.pi

ell=Ellipse((lupus3_table["H_est"],lupus3_table["sigma_est"]),ell_size[0],ell_size[1],ell_angle)
ell.set_alpha(0.3)


ax.add_artist(ell)
ax.text(lupus3_table["H_est"],lupus3_table["sigma_est"],1,ha="center",va="center",bbox=bbox_props)

ax.set_xticks(np.linspace(0.1,0.9,5))
ax.set_yticks(np.linspace(0.7,2.3,5))

ax.set_xlabel(r"$H$")
ax.set_ylabel(r"$\sigma$")

fig_4.tight_layout()
fig_4.savefig("H_sigma_plot.png",dpi=150)




"""
taurus_table=taurus.make_table_row()
chai_table=chai.make_table_row()
ic348_table=ic348.make_table_row()
lupus3_table=lupus3.make_table_row()
ophiuchus_table=ophiuchus.make_table_row()
ic2391_table=ic2391.make_table_row()
"""


"""

fig_1=plt.figure()

ax=fig_1.add_subplot(111)
ax.set_xlabel("sbar")
ax.set_ylabel("mbar")
ax.set

plot=ax.scatter(estimator_2d2d_2.test_data["s_mean"]/estimator_2d2d_2.test_data["R_hull"],estimator_2d2d_2.test_data["m_mean"]/estimator_2d2d_2.test_data["R_hull"],marker=".",
            c=estimator_2d2d_2.test_data["H"],cmap=plt.cm.inferno)



plt.figure(0)
plt.xlim([0.1,0.7])
plt.ylim([0.1,0.7])
plt.scatter(estimator_2d2d_2.test_data["s_mean"]/estimator_2d2d_2.test_data["R_hull"],estimator_2d2d_2.test_data["m_mean"]/estimator_2d2d_2.test_data["R_hull"],marker=".",
            c=estimator_2d2d_2.test_data["H"],cmap=plt.cm.inferno)
#plt.plot(np.linspace(0.,1.,100),np.zeros(100))
#plt.plot(np.linspace(0.,1.,100),np.linspace(0.,1.,100),linewidth=2.0,color="black")
plt.clim([0.,1.])
plt.colorbar(label="$H$")

plt.scatter(observations[:,0],observations[:,1],marker=".",s=300)
plt.xlabel("$sbar$")
plt.ylabel("$mbar$")
plt.show()

plt.figure(1)
plt.xlim([0.1,0.7])
plt.ylim([0.1,0.7])
plt.scatter(estimator_2d2d_2.test_data["s_mean"]/estimator_2d2d_2.test_data["R_hull"],estimator_2d2d_2.test_data["m_mean"]/estimator_2d2d_2.test_data["R_hull"],marker=".",
            c=estimator_2d2d_2.test_data["sigma"],cmap=plt.cm.inferno)
plt.scatter(observations[:,0],observations[:,1],marker=".",s=300)
#plt.plot(np.linspace(0.,1.,100),np.zeros(100))
#plt.plot(np.linspace(0.,1.,100),np.linspace(0.,1.,100),linewidth=2.0,color="black")
#plt.clim([0.5,2.5])
#plt.colorbar(label="$\sigma$")
#plt.xlabel("$H$")
#plt.ylabel("$H_{est}$")
plt.show()

plt.figure(2)
plt.xlim([0.1,0.7])
plt.ylim([0.1,0.7])
plt.scatter(fractal_cluster["s_mean"]/fractal_cluster["R_hull"],fractal_cluster["m_mean"]/fractal_cluster["R_hull"],marker=".",
            c=fractal_cluster["D"],cmap=plt.cm.inferno)
plt.scatter(radial_cluster["s_mean"]/radial_cluster["R_hull"],radial_cluster["m_mean"]/radial_cluster["R_hull"],marker=".",
            c=radial_cluster["alpha"],cmap=plt.cm.inferno)
plt.scatter(observations[:,0],observations[:,1],marker=".",s=300)
#plt.plot(np.linspace(0.,1.,100),np.zeros(100))
#plt.plot(np.linspace(0.,1.,100),np.linspace(0.,1.,100),linewidth=2.0,color="black")
#plt.clim([0.5,2.5])
#plt.colorbar(label="$\sigma$")
#plt.xlabel("$H$")
#plt.ylabel("$H_{est}$")
plt.show()





plt.figure(1)
plt.xlim([0.,1.])
plt.ylim([0.5,2.5])
plt.scatter(estimator_2d2d_2.test_data["H"],estimator_2d2d_2.test_data["sigma"],marker=".",cmap=plt.cm.inferno,
            c=estimator_2d2d_2.test_data["H_est"]-estimator_2d2d_2.test_data["H"])
plt.clim([-0.4,0.4])
plt.colorbar()
plt.show()


plt.figure(2)
plt.xlim([0.,1.])
plt.ylim([0.5,2.5])
plt.scatter(estimator_2d2d_2.test_data["H"],estimator_2d2d_2.test_data["sigma"],marker=".",cmap=plt.cm.inferno,
            c=estimator_2d2d_2.test_data["sigma_est"]-estimator_2d2d_2.test_data["sigma"])
plt.clim([-0.8,0.8])
plt.colorbar()
plt.show()
"""


"""


data_3d2d_2=pd.read_csv("data_3d2d_2.dat")
estimator_3d2d_2=param_estimator(data_3d2d_2)


estimator_3d2d_2.estimate_params(data_2d2d_2)
estimator_2d2d_2.estimate_params(data_3d2d_2)

plt.figure(1)
plt.scatter(data_2d2d_2["H"],data_2d2d_2["H_est"])
plt.show()


plt.figure(2)
plt.scatter(data_2d2d_2["sigma"],data_2d2d_2["sigma_est"])
plt.show()


"""

"""
plt.figure(1)
plt.xlim([0.,1.])
plt.ylim([0.,1.])
plt.scatter(estimator_2d2d_0.test_data["H_est"],estimator_2d2d_0.test_data["H"],marker=".",c=estimator_2d2d_0.test_data["sigma"],cmap=plt.cm.coolwarm)
#plt.plot(np.linspace(0.,1.,100),np.zeros(100))
plt.plot(np.linspace(0.,1.,100),np.linspace(0.,1.,100))
plt.show()

plt.figure(2)
plt.xlim([0.,1.])
plt.ylim([0.,1.])
plt.scatter(estimator_2d2d_1.test_data["H_est"],estimator_2d2d_1.test_data["H"],marker=".",c=estimator_2d2d_1.test_data["sigma"],cmap=plt.cm.coolwarm)
#plt.plot(np.linspace(0.,1.,100),np.zeros(100))
plt.plot(np.linspace(0.,1.,100),np.linspace(0.,1.,100))
plt.show()


plt.figure(3)
plt.xlim([0.,1.])
plt.ylim([0.,1.])
plt.scatter(estimator_2d2d_2.test_data["H"],estimator_2d2d_2.test_data["H_est"],marker=".",c=estimator_2d2d_2.test_data["sigma"],cmap=plt.cm.coolwarm)
#plt.plot(np.linspace(0.,1.,100),np.zeros(100))
plt.plot(np.linspace(0.,1.,100),np.linspace(0.,1.,100),linewidth=2.0,color="black")
plt.clim([0.5,2.5])
plt.colorbar(label="$\sigma$")
plt.xlabel("$H$")
plt.ylabel("$H_{est}$")
plt.show()
"""

"""
plt.figure(4)
plt.xlim([0.5,2.5])
plt.ylim([0.5,2.5])
plt.scatter(estimator_2d2d_0.test_data["sigma_est"],estimator_2d2d_0.test_data["sigma"],marker=".",c=estimator_2d2d_0.test_data["H"],cmap=plt.cm.coolwarm)
#plt.plot(np.linspace(0.5,2.5,100),np.zeros(100))
plt.plot(np.linspace(0.5,2.5,100),np.linspace(0.5,2.5,100))
plt.show()

plt.figure(5)
plt.xlim([0.5,2.5])
plt.ylim([0.5,2.5])
plt.scatter(estimator_2d2d_1.test_data["sigma_est"],estimator_2d2d_1.test_data["sigma"],marker=".",c=estimator_2d2d_1.test_data["H"],cmap=plt.cm.coolwarm)
#plt.plot(np.linspace(0.5,2.5,100),np.zeros(100))
plt.plot(np.linspace(0.5,2.5,100),np.linspace(0.5,2.5,100))


plt.figure(6)
plt.xlim([0.5,2.5])
plt.ylim([0.5,2.5])
plt.scatter(estimator_2d2d_2.test_data["sigma"],estimator_2d2d_2.test_data["sigma_est"],marker=".",c=estimator_2d2d_2.test_data["H"],cmap=plt.cm.coolwarm)
#plt.plot(np.linspace(0.5,2.5,100),np.zeros(100))
plt.plot(np.linspace(0.5,2.5,100),np.linspace(0.5,2.5,100),linewidth=2.0,color="black")
plt.clim([0,1.0])
plt.colorbar(label="$H$")
plt.xlabel("$\sigma$")
plt.ylabel("$\sigma_{est}$")
plt.show()
plt.show()
"""