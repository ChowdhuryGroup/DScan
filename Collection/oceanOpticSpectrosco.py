# -*- coding: utf-8 -*-

"""

UPDATED: 11-21-2019

DESCRIPTION:

    Read and write commands to SpectraPro HRS-300

    **NEW** Control Ocean Optics spectrometers 

    

"""

import serial as s

import seabreeze.spectrometers as sb #library for OceanOptics 

import time





class mono:   # create a monochromator class

    def comset(self,num):

        self.COMport=str(num)

        self.ser=s.Serial('COM'+str(self.COMport)) 

        self.ser.timeout=5

        self.ser.close()

    

    def oport(self): #check if COM port is open

        if self.ser.isOpen()==True:

            print('COM Port is already open.')

        else:

           self.ser.open()

        return;   

    

    def info(self): # check model and serial number of monochromator

      self.oport()

      self.ser.write('serial \r'.encode())

      monosn=self.ser.readline().decode()[:-4:]  #the [:-4:] removes the last 4 indices of the string, which includes the 'ok'. This is optimal value.

      print('SN: '+monosn) 

      self.ser.write('model \r'.encode())

      monomod=self.ser.readline().decode()[:-4:]

      print('Model: '+monomod)

      textout='Model: '+monomod+' || SN: '+monosn

      self.ser.close()

      return textout

    

    def gratings(self): #list installed gratings

      self.oport()

      self.ser.write('?gratings \r'.encode())

      monograt=self.ser.read_until('ok'.encode()).decode()[4:-3:]

      print('Installed Gratings: '+monograt)

      self.ser.close()

      return monograt

    

    def gratnum(self): # number for grating in use

      self.oport()

      self.ser.write('?grating \r'.encode())

      gratnum=self.ser.readline().decode()[1:-5:] 

      print('Grating: '+gratnum)

      self.ser.close()

      return gratnum

  

    def gratinfo(self): # info about grating in use

      self.oport()

      gratnum=self.gratnum()

      monograt=self.gratings()

      ind1=monograt.find(gratnum)

      ind2=monograt.find('\r',ind1)

      gratinfo='Grating: '+monograt[ind1:ind2:]

      self.ser.close()

      return gratinfo          

        

        

    def state(self): # check wavelength and grating being used

      self.oport()

      self.ser.write('?nm \r'.encode())

      wl='Wavelength: '+self.ser.readline().decode()[:-4:]

      gratinfo=self.gratinfo()

      state=wl+' || '+gratinfo

      print(state)

      self.ser.close()

      return state;

    

    def setwl(self,wl): # set wavelength

        self.oport()

        self.ser.write((str(wl)+' goto \r').encode())

        print(self.ser.readline().decode())

        self.ser.write('?nm \r'.encode())

        print('Now set to '+self.ser.readline().decode()[:-4:])

        self.ser.close()

        return;



    def setgr(self,gr):

        self.oport()

        self.ser.write((str(gr)+' grating \r').encode())

        time.sleep(15)

        print(self.ser.readline().decode())

        print('IN PLACE: '+self.gratinfo())

        self.ser.close()

        return

        

class ocean:

    

    def __init__(self,sernum):

       self.sernum=sernum;

       self.spec=sb.Spectrometer.from_serial_number(sernum);
       #self.spec.trigger_mode(0)
       #sb.seabreeze.pyseabreeze.SeaBreezeThermoElectricFeature.enable_tec(True)

        

    def setinttime(self,num):

        num=num*1000

        self.spec.integration_time_micros(num);

        return



    

    def getspec(self):

        spectrum=self.spec.spectrum()

        return spectrum;
    
    def close(self):
        self.spec.close()
    

              

        

          

