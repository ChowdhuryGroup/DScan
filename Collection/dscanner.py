# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 14:43:12 2021

@author: Conrad Kuz
"""

import oceanOpticSpectrosco as spectro
import newfocusStage
import numpy as np
import matplotlib.pyplot as plt
import time
from tkinter.filedialog import asksaveasfilename
import csv
from datetime import datetime
from tqdm import tqdm



specSerialNumber = 'HR4P0326'
spectrum = spectro.ocean(specSerialNumber)
try:
    spectrum.setinttime(100)
except:
    print("except")
    spectrum.setinttime(100)


comPort = 'COM14'
stage = newfocusStage.smc100(comPort)

if input("capture background and fundamental? y/n: ") == "y":
    input("press enter to capture background: ")
    background = spectrum.getspec()
    plt.ion()
    plt.figure(1)
    plt.plot(background[0],background[1])
    plt.show()


    save = "n"
    fundamental = []

    while save != "y":
        input("press enter to capture fundamental: ")
        fundamental = list(spectrum.getspec())
        plt.figure(2)
        plt.clf()
        plt.plot(fundamental[0],fundamental[1])
        plt.show()
        plt.pause(0.05)
        plt.waitforbuttonpress()
        save = input("Save? y/n: ")



input("begin dscan press enter: ")

minPosition = 18
maxPosition = 18.92
numberPositions = 127 #use (2^n)-1
spectrometerIntegrationTime = 100 #ms
stepSize = (maxPosition-minPosition)/numberPositions
positions = np.arange(minPosition,maxPosition+stepSize,stepSize)
spectrum.setinttime(spectrometerIntegrationTime)

wavelengths = []
intensities = []
pbar = tqdm(positions, desc="Scanning")
for pos in pbar:
    stage.move_absolute(pos)
    stage.wait_till_done()
    pbar.set_description(f"Pos(mm): {pos:.3f}")
    time.sleep(.01)
    data = spectrum.getspec()
    wavelengths.append(data[0])
    intensities.append(data[1])

plt.ioff()
plt.figure(3)
[wls, poss] = np.meshgrid(np.array(wavelengths[1]),np.array(positions))
intensities = np.array(intensities)
plt.pcolormesh(wls, poss, intensities)
plt.xlim((350,450))
plt.show()

# --- Extract ridge (max intensity at each position) ---
wavelengths_arr = np.array(wavelengths)
intensities_arr = np.array(intensities)

ridge_wavelengths = []
left_index = np.abs(wavelengths_arr[0] - 350).argmin()
right_index = np.abs(wavelengths_arr[0] - 450).argmin()
for i in range(len(intensities_arr)):
    idx_max = np.argmax(intensities_arr[i][left_index:right_index]) + left_index
    ridge_wavelengths.append(wavelengths_arr[i][idx_max])

ridge_wavelengths = np.array(ridge_wavelengths)
positions_arr = np.array(positions)
left_index = ridge_wavelengths.argmin()
right_index = ridge_wavelengths.argmax()
ridge_wavelengths = ridge_wavelengths[left_index:right_index]
positions_arr = positions_arr[left_index:right_index]
# --- Fit line (position vs wavelength) ---
coeffs = np.polyfit(positions_arr, ridge_wavelengths, 1)
slope = coeffs[0]
intercept = coeffs[1]

print(f"Slope (nm/mm): {slope:.6f}")

# --- Generate fitted line ---
fit_line = slope * positions_arr + intercept

# --- Plot result ---
plt.figure(4)
plt.plot(positions_arr, ridge_wavelengths, 'o', label="Extracted ridge")
plt.plot(positions_arr, fit_line, '-', label=f"Fit (slope={slope:.4f})")
plt.xlabel("Position (mm)")
plt.ylabel("Wavelength (nm)")
plt.title("Streak Linear Fit")
plt.legend()
plt.show()

print("Scan Done!")
spectrum.close()
stage.close()
timestamp = datetime.now().isoformat()
filename = asksaveasfilename()
with open(filename, mode='w',newline='') as out_file:
    tsv = csv.writer(out_file,delimiter='\t')
    tsv.writerow(["# Timestamp", timestamp])
    tsv.writerow(["# minPosition", minPosition])
    tsv.writerow(["# maxPosition", maxPosition])
    tsv.writerow(["# numberPositions", numberPositions])
    tsv.writerow(["# integration_time_ms", spectrometerIntegrationTime])
    tsv.writerow(np.concatenate([[0],fundamental[0]]))
    tsv.writerow(np.concatenate([[0],fundamental[1]]))
    for i in range(len(positions)):
        tsv.writerow(np.concatenate([[positions[i]],intensities[i]]))
with open(filename[:-4]+"_background"+filename[-4:], mode = 'w', newline ='') as out_file:
    tsv = csv.writer(out_file,delimiter = '\t')
    tsv.writerow(["# Timestamp", timestamp])
    tsv.writerow(background[0])
    tsv.writerow(background[1])
    