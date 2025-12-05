#!/usr/bin/env python
# coding: utf-8

# In[10]:


import warnings
warnings.simplefilter('ignore')

import numpy as np
import subprocess
import lifesim
from spectres import spectres
import os
import matplotlib.pyplot as plt
from astropy import units as u
from matplotlib.pyplot import figure, show
import pandas

#correct for directory changes by lifesim
os.chdir('/Users/users/cesario/BULK/Retrievals-on-forward')

def read_dat_file_with_headers(file_path):
    data = np.loadtxt(file_path, comments='#')
    temperature = data[:, 0]
    pressure = data[:, 1]
    return temperature, pressure

def limits_function(file_path):
    data = np.loadtxt(file_path, comments='#')
    return tuple(data[:, i] for i in range(data.shape[1]))

def noack_files(filename, figure):
    data = np.loadtxt(f"Figure{figure}/{filename}",delimiter=',')
    return tuple(data[:,i] for i in range(data.shape[0]))

figure = str(input("Figure 2 or 3: "))

redox_state = float(input("Redox state: (6, 4.9, 3.8, 2.7, 1.6, 0.5 (and negatives)): "))
surface_temp = float(input("Surface temperature: (100, 190, 236, 281,327, 372, 418, 463, 509, 554, 600): "))

star_type = input("Orbiting star type (S, L or T): ")
distance = float(input("Distance star-planet system (pc): "))

IW = np.loadtxt(f"Figure{figure}/Figure{figure}_buffer.csv", delimiter=',')  # x-axis (IW values)

if figure == '2':
    T_surf = np.loadtxt(f"Figure{figure}/Figure{figure}_T_surf.csv", delimiter=',')  # y-axis (T values)
    
pressure_overall = np.loadtxt(f"Figure{figure}/Figure{figure}_P_atm.csv", delimiter=',')

pressure_H2 = noack_files(f"Figure{figure}_P_H2.csv",figure)
pressure_H2O = noack_files(f"Figure{figure}_P_H2O.csv",figure)
pressure_NH3 = noack_files(f"Figure{figure}_P_NH3.csv",figure)
pressure_N2 = noack_files(f"Figure{figure}_P_N2.csv",figure)
pressure_CH4 = noack_files(f"Figure{figure}_P_CH4.csv",figure)
pressure_CO2 = noack_files(f"Figure{figure}_P_CO2.csv",figure)

#shape: num_T x num_IW x 6 molecules
pressures = np.stack([pressure_H2, pressure_H2O, pressure_NH3, 
                      pressure_N2, pressure_CH4, pressure_CO2], axis=-1)

#IW and Temperature to print out:
IW_print = redox_state
T = surface_temp

if figure == '2':
    for i, iw in enumerate(IW):
        for j, temp in enumerate(T_surf):
            fractions = pressures[j, i, :]
            total_pressure = np.sum(fractions)
            if round(iw, 1) == IW_print and round(temp, 0) == T:
                overall_pressure = total_pressure
                print("Overall Pressure: ", total_pressure)
                H2 = fractions[0] / total_pressure
                H2O = fractions[1] / total_pressure
                NH3 = fractions[2] / total_pressure
                N2 = fractions[3] / total_pressure
                CH4 = fractions[4] / total_pressure
                CO2 = fractions[5] / total_pressure
                
elif figure == '3':
    for i, iw in enumerate(IW):
            fractions = pressures[-1, i, :]
            total_pressure = np.sum(fractions)
            if round(iw, 1) == IW_print:
                overall_pressure = total_pressure
                print("Overall Pressure: ", total_pressure)
                H2 = fractions[0] / total_pressure
                H2O = fractions[1] / total_pressure
                NH3 = fractions[2] / total_pressure
                N2 = fractions[3] / total_pressure
                CH4 = fractions[4] / total_pressure
                CO2 = fractions[5] / total_pressure
else:
    print(f"Error: Figure values has to be '2' or '3'")
            
print(F"H2:{H2}")
print(F"H2O:{H2O}")
print(F"NH3:{NH3}")
print(F"N2:{N2}")
print(F"CH4:{CH4}")
print(F"CO2:{CO2}")

if star_type == "S":
    Tstar=5777
    Rstar=1
    Dplanet=1
elif star_type == "L":
    Tstar=3390
    Rstar=0.4
    Dplanet=0.13
elif star_type == "T":
    Tstar=2566
    Rstar=0.12
    Dplanet=0.02
else:
    raise "NameError: Unrecognized star type"
    
file_string = f"""H2={H2}
H2O={H2O}
NH3={NH3}
N2={N2}
CH4={CH4}
CO2={CO2}

chemistry=.false.

*star --> sun at 1AU; sistem at 10pc
Tstar={Tstar}
Rstar={Rstar}
Dplanet={Dplanet}
distance={distance}

TeffP=2.7

*Earth radius and mass
Rp=0.089d0
Mp=0.00315d0

nr=50

lmin=0.1d0
lmax=20d0

specres=200

pmin=1d-8
pmax={overall_pressure}
Pp={overall_pressure}

scattering=.true.
scattstar=.true.

run3D=.true.
night2day=1

ComputeT=.true.
maxiter=30

cia1:file="/Users/users/cesario/ARCiS/Data/CIA/CO2-CO2_combined.cia"
cia2:file="/Users/users/cesario/ARCiS/Data/CIA/N2-N2_ordered.cia"
cia3:file="/Users/users/cesario/ARCiS/Data/CIA/N2-H2O_2018.cia"
cia4:file="/Users/users/cesario/ARCiS/Data/CIA/CO2-CH4_2018.cia"
cia5:file="/Users/users/cesario/ARCiS/Data/CIA/CO2-H2_2018.cia"
cia6:file="/Users/users/cesario/ARCiS/Data/CIA/H2-CH4_eq_2011.cia"

"""
print(f"Forward-Files/{star_type}_IW{round(IW_print, 0):.0f}_{surface_temp:.0f}.in")
#save to new .in file
with open(f"Forward-Files/{star_type}_IW{round(IW_print, 0):.0f}_{surface_temp:.0f}.in", "w") as f:
    f.write(file_string)

def replace_e_with_d_in_first_6_lines(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    for i in range(min(6, len(lines))):
        lines[i] = lines[i].replace('e', 'd').replace('E', 'D')

    with open(file_path, 'w') as f:
        f.writelines(lines)
         
    return

#replace exponential e with d, to match ARCiS
replace_e_with_d_in_first_6_lines(f"Forward-Files/{star_type}_IW{round(IW_print, 0):.0f}_{surface_temp:.0f}.in")
         
#run ARCiS forward model
subprocess.run(f"module load ARCiS && module load cfitsio && export OMP_NUM_THREADS=8 && ARCiS Forward-Files/{star_type}_IW{round(IW_print, 0):.0f}_{surface_temp:.0f}.in -o Forward-Files/NO_{star_type}_{round(IW_print, 0):.0f}_{int(surface_temp)} && echo 'Forward Model Successfully Generated!'", shell=True)

#lifesim OBSERVATION
emis = np.loadtxt(f"/Users/users/cesario/BULK/Retrievals-on-forward/Forward-Files/NO_{star_type}_{round(IW_print, 0):.0f}_{int(surface_temp)}/emis",comments='#')

integration_time = float(input("LIFE Observation / time (hours): "))
lambdaa = emis[:,0] * u.micron
wavelength = lambdaa.to(u.m)
jansky_flux = emis[:,1] * u.Jansky

#wavelength in angstrom for the unit conversion
lambdaa = lambdaa.to(u.angstrom)

photon_flux = (jansky_flux.value * 1.51e3)/lambdaa.value * 10e13 * u.ph / u.m**3 / u.s

#integration time
intime = integration_time #hours
intime *= 3600

#star observing parameters
dist = distance #parsec

#STAR parameters

#TRAPPIST-1=0, AD LEONIS=1, SUN=2
if star_type == 'T':
    star = 0
elif star_type == 'L':
    star = 1
elif star_type == 'S':
    star = 2

stemp = [2566,3390,5780][star] #Kelvin
sradius = [0.12,0.4,1][star] #Solar radii
slat = 0.78 #ecliptic latitude in radians 45deg = 0.78; 90deg = 1.57
orbit = [0.02,0.13,1][star] #AU
zodi = 3 #zodis, examples: 1,10,100,1000

#Instrument set-up, choose observing scenario
bus = lifesim.Bus()
#observing scenarios: optimistic, baseline, pessimistic
bus.data.options.set_scenario('baseline')
#spectral resolution
specres = 100
bus.data.options.set_manual(spec_res = specres)
#wavelength coverage
bus.data.options.set_manual(wl_min = 6)
bus.data.options.set_manual(wl_max = 16)
baseline_planet = True #optimizes baseline for planet

instrument = lifesim.Instrument(name='inst')
bus.add_module(instrument)
transm = lifesim.TransmissionMap(name='transm')
bus.add_module(transm)
#include the noise sources
exozodi = lifesim.PhotonNoiseExozodi(name='exo')
bus.add_module(exozodi)
localzodi = lifesim.PhotonNoiseLocalzodi(name='local')
bus.add_module(localzodi)
star_leak = lifesim.PhotonNoiseStar(name='star')
bus.add_module(star_leak)
#connect the modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('star', 'transm'))

planet_flux = [wavelength, photon_flux]

SNR = 0
n = orbit * u.AU; n = n.to(u.pc) #radius of orbit
ang = (n.value/dist) * (360/(2*np.pi)) * 3600

while SNR<10:
    bins_snr, pflux, noyz = instrument.get_spectrum(temp_s=stemp, radius_s=sradius, lat_s=slat, z=zodi, distance_s=dist, angsep=ang, flux_planet_spectrum=planet_flux, integration_time=intime)
    SNR = np.sqrt(np.sum(bins_snr[1]**2))
    print(f"S/N: {SNR} \n Integration time: {intime/3600:.1f} hours \n")
    intime += 3600

print(f"S/N of observation: ",SNR,f" Integration time: {intime/3600:.1f} hours")

planet_flux_r = spectres(new_wavs=instrument.data.inst['wl_bin_edges'],spec_wavs=planet_flux[0].value,spec_fluxes=planet_flux[1].value,edge_mode=True)

#noise distribution
#calculate simulated noise
ratio = planet_flux_r / bins_snr[1]
noise = np.random.normal(0,ratio,size=pflux.shape)

bins = bins_snr[0] * u.m
bins = bins.to(u.micron)

width = (bins.value[-1] - bins.value[0])/specres * u.micron
resolution = bins.value/width.value

#first we have to convert the emission spectra back to Jansky
preconv = bins_snr[0] * u.meter
conversion_lambda = preconv.to(u.angstrom)

mask = (preconv.to(u.micron).value > 6) & (preconv.to(u.micron).value < 16)

dat_emission = ((planet_flux_r * conversion_lambda.value)/(1.5e3 * 10e13))
dat_emission = dat_emission[mask]
dat_wavelength = preconv.to(u.micron)[mask]
dat_error = (noise * conversion_lambda.value)/(1.5e3 * 10e13)
dat_error = dat_error[mask]
dat_specres = resolution[mask]

filename = f"Fluxes/NO_{star_type}_{round(IW_print, 0):.0f}_{int(surface_temp)}"

data = np.column_stack((dat_wavelength.value, dat_emission, dat_error, dat_specres))
np.savetxt(filename, data, fmt='%.6e', delimiter=' ')

subprocess.run(f"echo 'Flux file post-observation created at Fluxes/{filename}. Now go retrieve it!'", shell=True)

#start the retrieval

subprocess.run(f"tmux", shell=True)

#&& ARCiS retrieval_CIA_partial.in -o NO_{star_type}_{round(IW_print, 0):.0f}_{int(surface_temp)}_test -s obs1:file=Fluxes/NO_{star_type}_{round(IW_print, 0):.0f}_{int(surface_temp)}

# 
