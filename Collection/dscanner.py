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




specSerialNumber = 'HR4P0326'
spectrum = spectro.ocean(specSerialNumber)
try:
    spectrum.setinttime(100)
except:
    print("except")
    spectrum.setinttime(100)


comPort = 'COM1'
stage = newfocusStage.smc100(comPort)

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

minPosition = 12
maxPosition = 15
numberPositions = 127 #use (2^n)-1
stepSize = (maxPosition-minPosition)/numberPositions
positions = np.arange(minPosition,maxPosition+stepSize,stepSize)

wavelengths = []
intensities = []

for pos in positions:
    stage.move_absolute(pos)
    stage.wait_till_done()
    print(pos)
    time.sleep(.01)
    data = spectrum.getspec()
    wavelengths.append(data[0])
    intensities.append(data[1])

plt.figure(3)
[wls, poss] = np.meshgrid(np.array(wavelengths[1]),np.array(positions))
intensities = np.array(intensities)
plt.pcolormesh(wls, poss, intensities)
plt.show()
print("Scan Done!")
spectrum.close()
stage.close()

filename = asksaveasfilename()
with open(filename, mode='w',newline='') as out_file:
    tsv = csv.writer(out_file,delimiter='\t')
    tsv.writerow(np.concatenate([[0],fundamental[0]]))
    tsv.writerow(np.concatenate([[0],fundamental[1]]))
    for i in range(len(positions)):
        tsv.writerow(np.concatenate([[positions[i]],intensities[i]]))
with open(filename[:-4]+"_background"+filename[-4:], mode = 'w', newline ='') as out_file:
    tsv = csv.writer(out_file,delimiter = '\t')
    tsv.writerow(background[0])
    tsv.writerow(background[1])
    