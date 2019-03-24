import os
import time
import glob

def findMostRecentFile():
    #find most recent file in log folder
    fName=max(glob.iglob('*.dat'),key=os.path.getmtime)
    return fName

def main():
	filename = findMostRecentFile()
	f = open(filename)
	line = f.readlines()[-1].split(',')
	f.close()
	header_list = ["Filename","Time Stamp","Status","Temperature","Magnetic Field","Sample Position","Bridge 1 Resistance","Bridge 1 Excitation","Bridge 2 Resistance","Bridge 2 Excitation","Bridge 3 Resistance","Bridge 3 Excitation","Bridge 4 Resistance","Bridge 4 Excitation","Signal 1 Vin","Signal 2 Vin","Digital Inputs","Drive 1 Iout","Drive 1 Ipower","Drive 2 Iout","Drive 2 Ipower","Pressure","Br_ch#2","Map 21","Map 22","Map 23","Map 24","Map 25","Map 26","Map 27","Map 28","Map 29","System Status","Setpoint","Magnet","Valve Position","CIG/Cernox","CIG Excit","Plat Therm","Plat Excit","Neck Therm","Neck Excit","Spare Therm","Spare Excit","Flow Rate","Annulus Vacuum","Impedance","Block Htr Cur","Block Htr Pwr","Neck Htr Cur","Neck Htr Pwr","Annulus Pres","System Pres","Tableless Map53","Tableless Map54","Redirected 55","Redirected 56","Redirected 57","Redirected 58","Redirected 59","Redirected 60","Redirected 61","Seq Stat","Sys Temp","He Level","Map 67","Map 68","Map 69","low res (A) read","Hi res (A) read","Mag (V) read","Mag adj/sc(I)set","Sw heater state","DAC <I>","Therm v","KRPMs","Pump temp","Pump power","High Vacuum mode","Flapper state","Regen htr flag","iso flag","BlockHeater","NeckHeater","Last 2 eng st","Hires enable","Mag (V) set","Map 89","Map 90","Map 91","Map 92","Map 93","Sytem Status","EverCool Mode","Low Pressure1","Low Pressure2","High Pressure","Pressure PID On/Off","Packed Valve","Recirc Valve","Exhaust Valve","Heater On/Off Read","Heater Votls Read","Compressor State","Compressor on time","Sytem Status","HeaterPower","HeaterPower","Thermometer 1","Thermometer 2","therm config","CompressorTimeoutTimer","MAP 21","MAP 22","MAP 23","MAP 24","MAP 25","MAP 26","MAP 27","MAP 28","MAP 29","MAP 30"]
	outfile = {}
	outfile["Filename"] = filename
	for index in range(len(line)):
		try:
			outfile[header_list[index]] = float(line[index])
		except:
			outfile[header_list[index]] = line[index]
main()