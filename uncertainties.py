import numpy as np
from matplotlib import pyplot as plt
from clusterfrac.estimator import param_estimator
import pandas as pd
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText


# set up models

data_2d2d_0=pd.read_csv("data_2d2d_0.dat")
estimator_2d2d_0=param_estimator(data_2d2d_0)

data_2d2d_1=pd.read_csv("data_2d2d_1.dat")
estimator_2d2d_1=param_estimator(data_2d2d_1)

data_2d2d_2=pd.read_csv("data_2d2d_2.dat")
estimator_2d2d_2=param_estimator(data_2d2d_2)



fig_1=plt.figure(figsize=(12,6))
ax=[fig_1.add_subplot(2,3,i+1) for i in range(6)]


for i in range(3):
    ax[i].set_xlim(0,1)
    ax[i].set_ylim(0,1)
    
for i in range(3,6):
    ax[i].set_xlim(0.5,2.5)
    ax[i].set_ylim(0.5,2.5)



series_0=ax[0].scatter(estimator_2d2d_0.test_data["H"],estimator_2d2d_0.test_data["H_est"],marker=".",c=estimator_2d2d_0.test_data["sigma"],cmap=plt.cm.coolwarm)
ax[0].plot(np.linspace(0.,1.,10),np.linspace(0.,1.,10),color="black")
series_1=ax[1].scatter(estimator_2d2d_1.test_data["H"],estimator_2d2d_1.test_data["H_est"],marker=".",c=estimator_2d2d_1.test_data["sigma"],cmap=plt.cm.coolwarm)
ax[1].plot(np.linspace(0.,1.,10),np.linspace(0.,1.,10),color="black")
series_2=ax[2].scatter(estimator_2d2d_2.test_data["H"],estimator_2d2d_2.test_data["H_est"],marker=".",c=estimator_2d2d_2.test_data["sigma"],cmap=plt.cm.coolwarm)
ax[2].plot(np.linspace(0.,1.,10),np.linspace(0.,1.,10),color="black")

ax[0].set_yticks(np.linspace(0,1,6))
ax[0].set_xticks(np.linspace(0,0.8,5))
ax[0].set_ylabel(r"$H_\mathrm{est}$")
ax[0].set_xlabel(r"$H$")
ax[1].set_yticks([])
ax[1].set_xticks(np.linspace(0,0.8,5))
ax[1].set_xlabel(r"$H$")
ax[2].set_yticks([])
ax[2].set_xlabel(r"$H$")




at=AnchoredText(r"$32\leq N_\bigstar\leq99$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[0].add_artist(at)

rho_corr=np.corrcoef(estimator_2d2d_0.test_data["H_est"],estimator_2d2d_0.test_data["H"])[0,1]
rmse=np.sqrt(((estimator_2d2d_0.test_data["H_est"]-estimator_2d2d_0.test_data["H"])**2).mean())

print(rho_corr,rmse)

at=AnchoredText("RMSE = "+str(rmse.round(2))+"\nPCC = "+str(rho_corr.round(2)),prop=dict(size=10),frameon=True,loc=4)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[0].add_artist(at)



at=AnchoredText(r"$100\leq N_\bigstar\leq315$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[1].add_artist(at)

rho_corr=np.corrcoef(estimator_2d2d_1.test_data["H_est"],estimator_2d2d_1.test_data["H"])[0,1]
rmse=np.sqrt(((estimator_2d2d_1.test_data["H_est"]-estimator_2d2d_1.test_data["H"])**2).mean())

print(rho_corr,rmse)


at=AnchoredText("RMSE = "+str(rmse.round(2))+"\nPCC = "+str(rho_corr.round(2)),prop=dict(size=10),frameon=True,loc=4)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[1].add_artist(at)

at=AnchoredText(r"$316\leq N_\bigstar\leq999$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[2].add_artist(at)

rho_corr=np.corrcoef(estimator_2d2d_2.test_data["H_est"],estimator_2d2d_2.test_data["H"])[0,1]
rmse=np.sqrt(((estimator_2d2d_2.test_data["H_est"]-estimator_2d2d_2.test_data["H"])**2).mean())

print(rho_corr,rmse)

at=AnchoredText("RMSE = "+str(rmse.round(2))+"\nPCC = "+str(rho_corr.round(2)),prop=dict(size=10),frameon=True,loc=4)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[2].add_artist(at)


series_0.set_clim(0.5,2.5)
series_1.set_clim(0.5,2.5)
series_2.set_clim(0.5,2.5)

series_3=ax[3].scatter(estimator_2d2d_0.test_data["sigma"],estimator_2d2d_0.test_data["sigma_est"],marker=".",c=estimator_2d2d_0.test_data["H"],cmap=plt.cm.coolwarm)
ax[3].plot(np.linspace(0.5,2.5,10),np.linspace(0.5,2.5,10),color="black")
series_4=ax[4].scatter(estimator_2d2d_1.test_data["sigma"],estimator_2d2d_1.test_data["sigma_est"],marker=".",c=estimator_2d2d_1.test_data["H"],cmap=plt.cm.coolwarm)
ax[4].plot(np.linspace(0.5,2.5,10),np.linspace(0.5,2.5,10),color="black")
series_5=ax[5].scatter(estimator_2d2d_2.test_data["sigma"],estimator_2d2d_2.test_data["sigma_est"],marker=".",c=estimator_2d2d_2.test_data["H"],cmap=plt.cm.coolwarm)
ax[5].plot(np.linspace(0.5,2.5,10),np.linspace(0.5,2.5,10),color="black")


series_3.set_clim(0,1)
series_4.set_clim(0,1)
series_5.set_clim(0,1)

ax[3].set_yticks(np.linspace(0.5,2.5,5))
ax[3].set_xticks(np.linspace(0.5,2.0,4))
ax[3].set_ylabel(r"$\sigma_\mathrm{est}$")
ax[3].set_xlabel(r"$\sigma$")
ax[4].set_yticks([])
ax[4].set_xlabel(r"$\sigma$")
ax[4].set_xticks(np.linspace(0.5,2.0,4))
ax[5].set_yticks([])
ax[5].set_xlabel(r"$\sigma$")

at=AnchoredText(r"$32\leq N_\bigstar\leq99$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[3].add_artist(at)

rho_corr=np.corrcoef(estimator_2d2d_0.test_data["sigma_est"],estimator_2d2d_0.test_data["sigma"])[0,1]
rmse=np.sqrt(((estimator_2d2d_0.test_data["sigma_est"]-estimator_2d2d_0.test_data["sigma"])**2).mean())

print(rho_corr,rmse)


at=AnchoredText("RMSE = "+str(rmse.round(2))+"\nPCC = "+str(rho_corr.round(2)),prop=dict(size=10),frameon=True,loc=4)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[3].add_artist(at)

at=AnchoredText(r"$100\leq N_\bigstar\leq315$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[4].add_artist(at)

rho_corr=np.corrcoef(estimator_2d2d_1.test_data["sigma_est"],estimator_2d2d_1.test_data["sigma"])[0,1]
rmse=np.sqrt(((estimator_2d2d_1.test_data["sigma_est"]-estimator_2d2d_1.test_data["sigma"])**2).mean())

print(rho_corr,rmse)

at=AnchoredText("RMSE = "+str(rmse.round(2))+"\nPCC = "+str(rho_corr.round(2)),prop=dict(size=10),frameon=True,loc=4)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[4].add_artist(at)

at=AnchoredText(r"$316\leq N_\bigstar\leq999$",prop=dict(size=10),frameon=True,loc=2)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[5].add_artist(at)

rho_corr=np.corrcoef(estimator_2d2d_2.test_data["sigma_est"],estimator_2d2d_2.test_data["sigma"])[0,1]
rmse=np.sqrt(((estimator_2d2d_2.test_data["sigma_est"]-estimator_2d2d_2.test_data["sigma"])**2).mean())

print(rho_corr,rmse)

at=AnchoredText("RMSE = "+str(rmse.round(2))+"\nPCC = "+str(rho_corr.round(2)),prop=dict(size=10),frameon=True,loc=4)
at.patch.set_boxstyle("round, pad=0.,rounding_size=0.2")
ax[5].add_artist(at)

fig_1.tight_layout()
fig_1.subplots_adjust(wspace=0,right=0.88)



pos=ax[2].get_position()
cbar_ax = fig_1.add_axes([0.90, pos.y0, 0.02, pos.y1-pos.y0])
fig_1.colorbar(series_2, cax=cbar_ax,ticks=np.linspace(0.5,2.5,5),label=r"$\sigma$")

pos=ax[5].get_position()
cbar_ax = fig_1.add_axes([0.90, pos.y0, 0.02, pos.y1-pos.y0])
fig_1.colorbar(series_5, cax=cbar_ax,ticks=np.linspace(0,1,6),label=r"$H$")

fig_1.savefig("uncertainties.png",dpi=150)

