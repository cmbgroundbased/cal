import pycal
from mpi4py import MPI
import numpy as np
import pandas as pd
import yaml
import seaborn as sns
import matplotlib.pylab as plt
import datetime

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# Define a CES scanning strategy of the line of sight of the telescope.

year=2020
month=4
day=28
hour=20
primes=0
seconds=0

data=datetime.datetime(year, month, day, hour, primes, seconds)
timestamp=data.timestamp()

weather=pycal.Weather("weather_STRIP.fits")
weather.set_time(int(timestamp))

# H20 in atmpsphere
pwv=weather.pwv
ice=weather.ice_water
liq=weather.liquid_water
hum=weather.humidity

# Temperaure and pressure
t0=weather.surface_temperature
t10=weather.air_temperature
p0=weather.surface_pressure

# Wind speed and direction
wx=weather.west_wind
wy=weather.south_wind
w = np.sqrt(wx**2 + wy**2)
w_dir = np.arctan2(wy, wx)

# Sky patch parameters
azmin=(0/180)*np.pi
azmax=(359/180)*np.pi
elmin=(60/180)*np.pi
elmax=(80/180)*np.pi
tmin=0.0
tmax_sim=1000.0
tmax_tod=tmax_sim
lmin_center=0.01
lmin_sigma=0.001
lmax_center=10
lmax_sigma=10
w_center=w
w_sigma=0
wdir_center=w_dir
#wdir_center=np.pi/2
wdir_sigma=0
z0_center=2000
z0_sigma=0
T0_center=t10
T0_sigma=0
zatm=40000.0
zmax=2000.0
xstep=100.0
ystep=100.0
zstep=100.0
nelem_sim_max=100000
verbosity=1
comm=comm
key1=0
key2=2**32+1+1
counterval1=0
counterval2=0
cachedir="atm_cache_module_I"
rmin=0
rmax=10000


if rank == 0:
    env=pycal.Environment.get()
    print(env)

# Dummy scanning strategy
# samples frequency
fs_hz=20
# timestamp in seconds (now we simulate 1h of observations)
times=np.linspace(0,tmax_tod,int(tmax_tod)*fs_hz)

# Now we perform a really dummy CES at 72deg of elevation
# and with an AZ span of 10deg
el=np.ones(int(tmax_tod)*fs_hz)*(72/180)*np.pi
# Fixed telescope
#az=np.ones(int(tmax_tod)*fs_hz)*(270/180)*np.pi

# A dummy scanning strategy
azmin_tod=(60/180)*np.pi
azmax_tod=(315/180)*np.pi
az = (0.5*(azmax_tod - azmin_tod)*np.sin((2*np.pi/500)*times)) + (0.5*(azmax_tod + azmin_tod))

# import the focal plane unit
with open(r'/home/algebrato/.julia/dev/Stripeline/instrumentdb/strip_focal_plane.yaml') as file:
    focalplane=yaml.full_load(file)

detecotrs=focalplane['horns'].keys()
tod=dict()
pointings=dict()

theta_i = 0
theta_f = 90
bor = (((70/2)*np.sin(2*np.pi/150*times)+(70/2))/180)*np.pi

for i in detecotrs:
    directions=focalplane['horns'][i]['orientation']
    l=np.arctan(directions[0]/directions[2])
    u=np.arctan(directions[1]/directions[2])
    tod[i]={'signal':np.zeros(int(tmax_tod)*fs_hz), 'l':l, 'u':u}
    pointings[i]={'az':az+(l*np.cos(bor)-u*np.sin(bor)), 'el':el+(u*np.cos(bor)+l*np.sin(bor))}
    
obs=dict()
obs['name'] = "Strip_1H_Full_Focal_Plane"
obs['tod'] = tod
obs['pointings'] = pointings

data=pycal.Data(comm)

data.obs.append(obs)

# Simulate the atmosphere evolution
atm=pycal.AtmSimMPI(azmin, 
                    azmax, 
                    elmin, 
                    elmax, 
                    tmin, 
                    tmax_sim, 
                    lmin_center, 
                    lmin_sigma, 
                    lmax_center, 
                    lmax_sigma, 
                    w_center, 
                    w_sigma, 
                    wdir_center, 
                    wdir_sigma, 
                    z0_center, 
                    z0_sigma, 
                    T0_center, 
                    T0_sigma, 
                    zatm, 
                    zmax, 
                    xstep, 
                    ystep, 
                    zstep, 
                    nelem_sim_max, 
                    verbosity, 
                    comm, 
                    key1, 
                    key2, 
                    counterval1, 
                    counterval2, 
                    cachedir, 
                    rmin, 
                    rmax
)

err = atm.simulate(True)
comm.Barrier()

# for i in detecotrs:
#     atm.observe(times, data.obs[0]['pointings'][i]['az'], data.obs[0]['pointings'][i]['el'], data.obs[0]['tod'][i]['signal'], -1)
#     comm.Barrier()

# det_ob=data.obs[0]['tod']
# tod_matrix=dict()

# for i in detecotrs:
#     tod_matrix[i]=det_ob=data.obs[0]['tod'][i]['signal']

# matrix = pd.DataFrame(tod_matrix)

# matrix.to_csv(r"TOD_matrix.csv",index = False, header=True)


# MTX_CX=matrix.corr()

# plt.figure(1, figsize=(20,8.7))
# ax=sns.heatmap(MTX_CX, annot=False, cmap='coolwarm', vmin=.5, vmax=1, center=0.75, square=True)
# ax.set(xlabel='Detectors name', ylabel='Detectors name', title="Correlation Level")
# plt.show()

