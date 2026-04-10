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

minPosition = 17.5
maxPosition = 19.4
numberPositions = 63 #use (2^n)-1
spectrometerIntegrationTime = 150 #ms
spectrums_to_average = 1

stage = newfocusStage.smc100(comPort)
spectrum.setinttime(spectrometerIntegrationTime)

if input("capture background and fundamental? y/n: ") == "y":
    input("press enter to capture background: ")
    background = np.empty_like(spectrum.getspec())
    for _ in range(spectrums_to_average):
        background += spectrum.getspec()
    background = background / spectrums_to_average
    plt.ion()
    plt.figure(1)
    plt.plot(background[0],background[1])
    plt.show()


    save = "n"
    fundamental = np.empty_like(spectrum.getspec())

    while save != "y":
        input("press enter to capture fundamental: ")
        for _ in range(spectrums_to_average):
             fundamental += spectrum.getspec()
             #time.sleep(spectrometerIntegrationTime/1000)
        fundamental = fundamental / spectrums_to_average
        plt.figure(2)
        plt.clf()
        plt.plot(fundamental[0],fundamental[1])
        plt.show()
        plt.pause(0.05)
        plt.waitforbuttonpress()
        save = input("Save? y/n: ")



input("begin dscan press enter: ")


stepSize = (maxPosition-minPosition)/numberPositions
positions = np.arange(minPosition,maxPosition+stepSize,stepSize)


wavelengths = []
intensities = []
pbar = tqdm(positions, desc="Scanning")
for pos in pbar:
    stage.move_absolute(pos)
    stage.wait_till_done()
    pbar.set_description(f"Pos(mm): {pos:.3f}")
    time.sleep(.01)
    total_spectrum = np.zeros_like(spectrum.getspec())

    for _ in range(spectrums_to_average):
        total_spectrum += spectrum.getspec()
    total_spectrum = total_spectrum / spectrums_to_average
    wavelengths.append(total_spectrum[0])
    intensities.append(total_spectrum[1])

plt.ioff()
plt.figure(3)
[wls, poss] = np.meshgrid(np.array(wavelengths[1]),np.array(positions))
intensities = np.array(intensities)
plt.pcolormesh(wls, poss, intensities)
plt.xlim((350,450))
plt.show()

while True:
    # --- Show plot and collect clicks ---
    plt.figure(5)
    plt.pcolormesh(wls, poss, intensities)
    plt.xlim((350, 450))
    plt.title("Click TWO points to define slope")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Position (mm)")

    points = plt.ginput(2, timeout=-1)
    plt.close()

    # Make sure user actually clicked twice
    if len(points) < 2:
        print("You didn't select two points. Try again.")
        continue

    (x1, y1), (x2, y2) = points

    # Avoid divide-by-zero
    if y2 == y1:
        print("Points have same position (vertical line). Try again.")
        continue

    slope = (y2-y1)/(x2-x1)

    print(f"\nSelected points:")
    print(f"  Point 1: (λ={x1:.3f} nm, pos={y1:.3f} mm)")
    print(f"  Point 2: (λ={x2:.3f} nm, pos={y2:.3f} mm)")
    print(f"Slope (mm/nm): {slope:.6f}")

    # --- Show preview ---
    plt.figure(6)
    plt.pcolormesh(wls, poss, intensities)
    plt.scatter([x1, x2], [y1, y2], color='red')

    plt.plot([x1, x2], [y1, y2], 'r--', label=f"Slope={slope:.4f}")
    plt.xlim((350, 450))
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Position (mm)")
    plt.legend()
    plt.title("Preview (close window to continue)")
    plt.show()

    # --- Ask user if they want to keep it ---
    user_input = input("Accept this slope? (y/n): ").strip().lower()
    if user_input == 'y':
        break

print(f"\nFinal slope (nm/mm): {slope:.6f}")

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
    