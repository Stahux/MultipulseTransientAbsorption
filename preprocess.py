# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 17:19:09 2024

@author: Stanisław Niziński
"""

import sif_parser as sif
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import os
import re
import copy as cp

######################################################################
################THIS IS 1st STEP INITIAL PROCESSING CODE##############
################ it is meant to reduce *.sif files from ##############
################ the Andor camera and turn them into    ##############
################ managable smaller data files used later##############
################ for the data analysis. It is not meant ##############
################ to be fast or elegant. It is meant to  ##############
################ do its job correctly. It is always work##############
################ in progress and tailored for this      ##############
################ specific experiment.                   ##############
######################################################################

def newWavelengthGrid(old_grid, new_grid, old_surface):
    #project to new wavelength grid, be careful to make nice grid
    #each new point must cover a few old points
    #and it cannot span points outside old grid!

    new_surface = np.zeros([old_surface.shape[0],new_grid.shape[0]])
    counters = np.zeros(new_grid.shape[0])
    delta = new_grid[1] - new_grid[0]
    min_ind = np.argmin(np.abs(old_grid-new_grid[0]+delta/2))
    max_ind = np.argmin(np.abs(old_grid-new_grid[-1]-delta/2))
    
    for i in range(min_ind,max_ind+1): #find closest old kinetcs
        k = np.argmin(np.abs(old_grid[i]-new_grid))
        new_surface[:,k] += old_surface[:,i]
        counters[k] += 1
        
    if(any(counters == 0)):
        print(counters)
        raise Exception("Trying to project on bad grid! Errror!")
    
    for i in range(new_grid.shape[0]):
        new_surface[:,i] /= counters[i]
            
    return new_surface


################### LOAD DATA ###################
#################################################
directory = '10-07-2024 OCP 74' #here specify dataset that should be processed
#################################################
#################################################
#################################################

#save drop this code to out directory to to know later what was used to process this data
file_code = open("preprocess.py", "r")
file_codeout = open(directory+"\\preprocess_copy.py", "w")
code_text = file_code.read()
file_codeout.write(code_text)
file_code.close()
file_codeout.close()

#load merged data arrays (all data with empty regions cut and wavelengths bucket-averaged)
#or create merged data arrays (if they are not there yet...)
try:
    
    print("Trying to load sorted arrays...")
    
    all_kinetics = np.load(directory+"\\all_kinetics.npy")   
    numbers = np.load(directory+"\\numbers.npy")
    delays = np.load(directory+"\\delays.npy")
    wavelengths = np.load(directory+"\\wavelengths.npy")
    #wavelengths_notch = np.load(directory+"\\wavelengths_notch.npy")
    periods = np.load(directory+"\\periods.npy")
    periodpattern = np.tile(periods, 90)
    
    print("Loaded sorted arrays!")
    
except:
    
    print("Failed to load sorted arrays, so trying to load from sif files!")

    periods = np.array([0.00300, #these are t_delays
                        0.01000, #t_delays must be the same like in experiment!
                        0.03000,
                        0.10000,
                        0.30000,
                        1.00000,
                        3.00000,
                        10.0000,
                        30.0000,
                        100.000])
    periodpattern = np.tile(periods, 90)
    
    all_kinetics = [] #drop all data here
    timestamps = [] #in seconds, to check if spacing makes sense
    numbers = [] #numbers of the sif files. check ordering, sort eventually
    
    #calc it based on first data file
    timedelta = None #temporal distance between two timepoints
    wavelengths_detector = None #grid as registered by detector - default setting of Andor callibrated by vapor light
    delays = None
    #wavelengths = np.linspace(410,778,int((778-410)/2+1)) #grid after projection to new grid (valid for MCWHF2)    
    #wavelengths = np.linspace(520,778,int((778-520)/2+1)) #grid after projection to new grid (valid for MCWHF2 + longpass)    
    wavelengths = np.linspace(470,778,int((778-470)/2+1)) #grid after projection to new grid (valid for MBB1F1 and no longpass)    
    
    lowest_number = np.inf
    
    #loop to load all files to the memory - very slow and not efficient for large datasets
    for filename in os.listdir(directory):
        if filename.endswith(".sif"):
            file_path = os.path.join(directory, filename)
            
            file_str = re.findall(r'(\d+) (\d{4})_(\d{2})\.(\d{2})\.(\d{2})_(\d+).sif', filename)
            print(file_str)
            
            day = int(file_str[0][0])
            hour = int(file_str[0][0+2])
            minute = int(file_str[0][1+2])
            second = int(file_str[0][2+2])
            number = int(file_str[0][3+2])
            
            #if(number > 1900): #uncomment to skip some
            #    continue
            
            if(number < lowest_number):
                lowest_number = number
            
            second_count = second + 60*minute + 60*60*hour + 60*60*24*day
            
            #if(number > 400): #uncomment to skip some
            #    continue
            
            file = sif.np_open(file_path)
        
            
            if(timedelta == None): #do it only for the first file!
                timedelta = (file[1]["timestamp_of_1"]-file[1]["timestamp_of_0"])/1000000
                delays = timedelta*np.arange(1,file[0][:,0,0].shape[0]+1)
                wavelengths_detector = np.arange(1,file[0].shape[2]+1)*0.34272949+232.342048
                #above callibrated based on HgAr lamp spectral lines
  
            #project to less dense wavelength grid
            projected_surface = newWavelengthGrid(wavelengths_detector, wavelengths, file[0][:,0,:] )

            all_kinetics.append(projected_surface)
            timestamps.append(second_count)
            numbers.append(number)

    
    ################### SORT DATA ###################
    
    numbers = np.array(numbers)
    timestamps = np.array(timestamps)
    
    #sort this stuff
    ordering = np.argsort(numbers)
    numbers = numbers[ordering]
    timestamps = timestamps[ordering]
    
    all_kinetics_new = [0]*len(all_kinetics) #overcomplicated, but works for lists too
    for i in range(len(all_kinetics)):
        all_kinetics_new[i] = all_kinetics[ordering[i]]
    all_kinetics = all_kinetics_new
    
    if(lowest_number != 1): #fix for problem when first number is > 1
        print("First number of spectr. is %i, so I will subtract to have 1 here." % lowest_number)
        numbers = 1 + numbers - lowest_number
    
    
    #was needed before I added sorting, but can still be there
    if(not(np.all(numbers[1:]-numbers[:-1] >= 1))):
        raise Exception("Check data files ordering! Seems to be wrong!")
        
        
    all_kinetics = np.array(all_kinetics)
        
    ################### SAVE DATA ###################    
    
    np.save(directory+"\\all_kinetics", all_kinetics)   
    np.save(directory+"\\numbers", numbers)
    np.save(directory+"\\delays", delays)
    np.save(directory+"\\wavelengths", wavelengths)
    np.save(directory+"\\periods", periods)
    
    if(False): #mark notch filter presence?
        n_min = np.argmin(np.abs(wavelengths-508))
        n_max = np.argmin(np.abs(wavelengths-518))
        wavelengths_notch = cp.deepcopy(wavelengths)
        wavelengths_notch[n_min:n_max+1] = np.nan    
        np.save(directory+"\\wavelengths_notch", wavelengths_notch)
        
        
################### DIAGNOSTICS SEQUENCE BELOW ###################           
        
#firsly check if darks are in right places
#so one will know is there are problems with experiment sequence
isdark = np.full(numbers.shape, True, dtype=bool)
for i in range(isdark.shape[0]):
    avg_tmp = np.average(all_kinetics[i,:,:])
    if(avg_tmp > 2000): #if above threshold, then it is probably not dark
        isdark[i] = False
        if(numbers[i] % 4 == 0): #complain if light is where dark is expected
            raise Exception("Kinetic with number %i has avg signal of %0.1f but it should be dark!" % (numbers[i],avg_tmp))        
            
    else: #below threshold, so dark
        if(numbers[i] % 4 != 0): #complain if there is no light but there should be probe
            raise Exception("Kinetic with number %i has avg signal of %0.1f but shouldn't be dark!" % (numbers[i],avg_tmp))        
            

#find places where laser hit the sample (signal spikes) and determine if they are in right places...
#so this function is just to see 2 or more adjacent spectra with laser scattering as one
def reduce_true_islands_1d(arr): #helper func
    arr = np.array(arr, dtype=bool)
    reduced_arr = np.zeros_like(arr, dtype=bool)
    
    prev_true = False
    for i in range(len(arr)):
        if arr[i] and not prev_true:
            reduced_arr[i] = True
            prev_true = True
        elif not arr[i]:
            prev_true = False
    
    return reduced_arr

k_laser = np.argmin(np.abs(wavelengths-512)) #some point where laser scattering can be found
k_sref = np.argmin(np.abs(wavelengths-750)) #where things should be relatively stable

#below code just checks if sequence looks correct based on laser scattering times
#it just warns me in case in case something looks suspiciuos, then i can take a closer look
for i in range(numbers.shape[0]):
    tmp = all_kinetics[i,:,k_laser]/all_kinetics[i,:,k_sref]
    avg_tmp = np.average(all_kinetics[i,:,k_laser]/all_kinetics[i,:,k_sref])
    #MAY BE NECESSARY TO MODIFY THRESHOLD BELOW - 1.2 currently
    laser_present = tmp/avg_tmp > 1.2 #true for places with exceptionally high light -> laser scattering
    laser_true = reduce_true_islands_1d(laser_present)
    
    bin_no = int(np.floor((numbers[i]-1)/4) % periods.shape[0])
    pulse2d = 2.490 #xpected 2nd pulse
    pulse1d = pulse2d - periods[bin_no] #xpected 1st pulse
    
    laser_delays = []
    for bv in range(len(laser_true)): #determine laser times
        if(laser_true[bv]):
            laser_delays.append(delays[bv])
            
    if(numbers[i] % 4 == 1): #there we should have only probe no laser pulses
        if(len(laser_delays) > 0):
            print("In kinetic no "+str(numbers[i])+" we should have no laser but i found some: ", laser_delays)
    elif(numbers[i] % 4 == 2): #there we should have one laser pulse close to delay 2.5s
        if(len(laser_delays) != 1):
            print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
        elif(np.abs(laser_delays[0]-pulse2d) > 0.03):
            print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
    elif(numbers[i] % 4 == 3): #there we should have max two laser pulses, one at 2.5s and optionally one before
        
        if(pulse1d <= 0.0): #only one pulse should be visible
            if(len(laser_delays) != 1):
                print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
                print("->It is strange because there should be only 2nd pulse visible, delay is ", periods[bin_no])  
            elif(np.abs(laser_delays[0]-pulse2d) > 0.03):
                print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
                print("->It is strange because there should be only 2nd pulse visible, delay is ", periods[bin_no])
        elif(periods[bin_no] < 0.01):
            if(len(laser_delays) != 1 and len(laser_delays) != 2):
                print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
                print("->It is strange because there should be 1 or closely spaced 2 pulses visible, delay is ", periods[bin_no])  
            elif(np.abs(laser_delays[0]-pulse1d) > 0.03):
                print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
                print("->It is strange because there should be at least one pulse visible, delay is: ", periods[bin_no])            
        else:
            if(len(laser_delays) != 2):
                print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
                print("->It is strange because there should be 2 pulses visible, delay is ", periods[bin_no])  
            elif(np.abs(laser_delays[1]-pulse2d) > 0.03 or np.abs(laser_delays[0]-pulse1d) > 0.03):
                print("In kinetic no "+str(numbers[i])+" I found strange laser pulses: ", laser_delays)
                print("->It is strange because there should be 2 pulses visible, delay is ", periods[bin_no])           
                
        


########################################################
################### PROCESSING BELOW ###################   
########################################################



################### AVG DATA + MEDIAN FILTER ###################

#function below is used to reject outlier cycles if they have unusually different 
#(usually low) signal. then histograms are plotted to ensure that everything
#is all right and proper point were rejected, and that dataset is still balanced
#(meaning that outliers are evenly scattered over t_delays and not too much data is rejected)
#Note that for some samples, for example crystalline ones this does not work very well. Works well for solutions.
def medianAbs(all_kinetics, numbers, periods = None, period_no = None, pulses_no = None):

    ddist = 3 #timepoint to look at
    itolerance = 0.05 #how far from median points to tolerate, early 0.02 was found to be good for solution data
        
    dshape = all_kinetics[0].shape
    
    kvect2_avg = np.zeros([dshape[0],dshape[1]])
    kvect1_avg = np.zeros([dshape[0],dshape[1]])
    kvect4_avg = np.zeros([dshape[0],dshape[1]])
    
    kvect2_avg_n = 0 #this is for 1hv/2hv - with pump
    kvect1_avg_n = 0 #this is for reference (only probe)
    kvect4_avg_n = 0 #this is for dark (no probe no pump)
    
    
    intdist1 = []
    intdist2 = []
    wdist = np.argmin(np.abs(700-wavelengths)) #intensity at 650 or 700 (large intensity still, no ocp band)

    #calc medians
    #if single pulse data
    if(pulses_no is None):
        for i in range(len(all_kinetics)):
            if(numbers[i] % 3 == 2): #1hv % 4 == 2
                intdist2.append(all_kinetics[i][ddist,wdist]) 
            elif(numbers[i] % 3 == 1): #0hv % 4 == 1
                intdist1.append(all_kinetics[i][ddist,wdist])
                
    #if two pulse data            
    else: 
        for i in range(len(all_kinetics)):
            bin_no = int(np.floor((numbers[i]-1)/4) % periods.shape[0])
            if(bin_no == period_no): #only take selected ones
                if(pulses_no == 1): #do kinetic for one pulse
                    if(numbers[i] % 4 == 2): #1hv % 4 == 2
                        intdist2.append(all_kinetics[i][ddist,wdist]) 
                elif(pulses_no == 2): #do kinetic for two pulses
                    if(numbers[i] % 4 == 3): #2hv % 4 == 3
                        intdist2.append(all_kinetics[i][ddist,wdist]) 
                    
                if(numbers[i] % 4 == 1): #0hv % 4 == 1
                    intdist1.append(all_kinetics[i][ddist,wdist])                
                
    
    median1 = np.median(intdist1)
    median2 = np.median(intdist2)
    
    intdist1_removed = [] #just to mark removed points
    intdist2_removed = []
    
    
    #if single pulse data
    if(pulses_no is None):
        for i in range(len(all_kinetics)):
            refint = all_kinetics[i][ddist,wdist] #to compare with median
            
            if(numbers[i] % 3 == 2): #1hv % 4 == 2
                if(np.abs(median2-refint) < itolerance*refint):
                    kvect2_avg[:,:] += all_kinetics[i][:,:]
                    kvect2_avg_n += 1
                    intdist2_removed.append(False)
                else:
                    intdist2_removed.append(True)
            elif(numbers[i] % 3 == 1): #0hv % 4 == 1
                if(np.abs(median1-refint) < itolerance*refint):
                    kvect1_avg += all_kinetics[i][:,:]
                    kvect1_avg_n += 1
                    intdist1_removed.append(False)
                else:
                    intdist1_removed.append(True)
            elif(numbers[i] % 3 == 0): #dark no probe no pump
                kvect4_avg += all_kinetics[i][:,:]
                kvect4_avg_n += 1
                
    #if two pulse data
    else: 
        for i in range(len(all_kinetics)):
            bin_no = int(np.floor((numbers[i]-1)/4) % periods.shape[0])
            refint = all_kinetics[i][ddist,wdist] #to compare with median
            
            if(bin_no == period_no): #only take selected ones
                if(pulses_no == 1): #do kinetic for one pulse
                    if(numbers[i] % 4 == 2): #1hv % 4 == 2
                        if(np.abs(median2-refint) < itolerance*refint):
                            kvect2_avg[:,:] += all_kinetics[i][:,:]
                            kvect2_avg_n += 1
                            intdist2_removed.append(False)
                        else:
                            intdist2_removed.append(True)
                            
                elif(pulses_no == 2): #do kinetic for two pulses
                    if(numbers[i] % 4 == 3): #2hv % 4 == 3
                        if(np.abs(median2-refint) < itolerance*refint):
                            kvect2_avg[:,:] += all_kinetics[i][:,:]
                            kvect2_avg_n += 1
                            intdist2_removed.append(False)
                        else:
                            intdist2_removed.append(True)
                    
                if(numbers[i] % 4 == 1): #0hv % 4 == 1
                    if(np.abs(median1-refint) < itolerance*refint):
                        kvect1_avg += all_kinetics[i][:,:]
                        kvect1_avg_n += 1 
                        intdist1_removed.append(False)
                    else:
                        intdist1_removed.append(True)
                
            if(numbers[i] % 4 == 0): #dark no probe no pump - just use all darks
                kvect4_avg += all_kinetics[i][:,:]
                kvect4_avg_n += 1
                
                
    if(True): #show histogram
        w = 5
        plt.figure(dpi=300)
        plt.title("Delay between pulses %.3fs and %i pulses used" % (periods[period_no],pulses_no))
        plt.hist(intdist1, bins=np.arange(median1*0.95, median1*1.05 + w, w), color=[1,0,0,0.5], label="0hv")
        plt.hist(intdist2, bins=np.arange(median1*0.95, median1*1.05 + w, w), color=[0,1,0,0.5], label="1hv/2hv")
        plt.legend()
        plt.xlabel("Intensity bins")
        plt.ylabel("Frequency of given intensity")
        plt.savefig(directory+"\\histogram_%i.png" % period_no)
        plt.show()
        
        intdist1_rm = np.copy(intdist1)
        for z in range(len(intdist1_removed)):
            if(intdist1_removed[z] == False):
                intdist1_rm[z] = np.nan
        intdist2_rm = np.copy(intdist2)
        for z in range(len(intdist2_removed)):
            if(intdist2_removed[z] == False):
                intdist2_rm[z] = np.nan      
        
        plt.figure(dpi=300)
        plt.title("Delay between pulses %.3fs and %i pulses used" % (periods[period_no],pulses_no))
        plt.plot(np.arange(1,len(intdist1)+1,1), intdist1, "r-", label="0hv")
        plt.plot(np.arange(1,len(intdist1)+1,1), intdist1_rm, "rx", label="0hv removed")
        plt.plot(np.arange(1,len(intdist2)+1,1), intdist2, "g-", label="1hv/2hv")
        plt.plot(np.arange(1,len(intdist2)+1,1), intdist2_rm, "gx", label="1hv/2hv removed")
        plt.legend()
        plt.xlabel("Repetition number")
        plt.ylabel("Intensity")
        plt.savefig(directory+"\\scanprogress_%i.png" % period_no)
        plt.show()
        
    
    print("Medians (ref, signal): ", median1, median2)
                    
                
    kvect2_avg = kvect2_avg/kvect2_avg_n
    kvect1_avg = kvect1_avg/kvect1_avg_n
    kvect4_avg = kvect4_avg/kvect4_avg_n
    
    
    Abs = np.zeros([dshape[0],dshape[1]])
    
    #finally calculate transient absorption
    Abs = -np.log10((kvect2_avg-kvect4_avg)/(kvect1_avg-kvect4_avg))
    
    print("Kinetics used (ref, signal, dark): ", kvect1_avg_n, kvect2_avg_n, kvect4_avg_n)
    print("Number of files is %i" % len(all_kinetics))
    
    return Abs




dshape = all_kinetics[0].shape

Abs1hv = np.zeros([periods.shape[0],dshape[0],dshape[1]])
Abs2hv = np.zeros([periods.shape[0],dshape[0],dshape[1]])

#call median filtering func and calculate transient absorptions
for bin_no in range(periods.shape[0]):
    Abs1hv[bin_no] = medianAbs(all_kinetics, numbers, periods = periods, period_no = bin_no, pulses_no = 1)
    Abs2hv[bin_no] = medianAbs(all_kinetics, numbers, periods = periods, period_no = bin_no, pulses_no = 2)



################### SAVE DATA ###################

#drop transient absorption data
np.save(directory+"\\Abs1hv", Abs1hv) 
np.save(directory+"\\Abs2hv", Abs2hv) 

################### PLOT DATA ###################

#uncomment to cut first second
t_cut = np.argmin(np.abs(delays-1.0))

#delays = delays[t_cut:]
#Abs1hv = Abs1hv[:,t_cut:,:]
#Abs2hv = Abs2hv[:,t_cut:,:]
dshape = Abs1hv[0].shape
#delays = delays-2.5
delays = delays-0.5


if(True): #subtract absorbance before the laser pulse
    for bin_no in range(periods.shape[0]):
        Abs1hv_zero = np.average(Abs1hv[bin_no,3:13,:], axis=0)
        Abs2hv_zero = np.average(Abs2hv[bin_no,3:13,:], axis=0)
        Abs1hv[bin_no,:,:] -= np.tile(Abs1hv_zero,(dshape[0],1))
        Abs2hv[bin_no,:,:] -= np.tile(Abs2hv_zero,(dshape[0],1))


#some important wavelengths
i_550 = np.argmin(np.abs(wavelengths-550))
i_500 = np.argmin(np.abs(wavelengths-500))
i_472 = np.argmin(np.abs(wavelengths-472))
i_438 = np.argmin(np.abs(wavelengths-438))
i_678 = np.argmin(np.abs(wavelengths-678))
i_520 = np.argmin(np.abs(wavelengths-520))


#plot data
for bin_no in range(periods.shape[0]):
     
    plt.figure(dpi=300)
    plt.title(directory.replace("OCP", "run no")+", delay "+str(periods[bin_no])+"s")
    plt.plot(delays, Abs1hv[bin_no,:,i_550]-0.0, label="1hv")
    plt.plot(delays, Abs2hv[bin_no,:,i_550]-0.0, label="2hv")
    #plt.ylim(-0.002,0.005)
    plt.grid(linestyle="dashed")
    plt.legend(loc="upper right")
    plt.xlabel("Time (s)")
    plt.ylabel("\u0394A (550nm)")
    plt.show()
    
    
    plt.figure(dpi=300)
    plt.title(directory.replace("OCP", "run no")+", delay "+str(periods[bin_no])+"s")
    plt.plot(wavelengths, Abs1hv[bin_no,579,:]-0.0, "b-", label="1hv %.3fs" % delays[579])
    plt.plot(wavelengths, Abs2hv[bin_no,579,:]-0.0, "b--", label="2hv %.3fs" % delays[579])
    plt.plot(wavelengths, Abs1hv[bin_no,767*2,:]-0.0, "g-", label="1hv %.3fs" % delays[767*2])
    plt.plot(wavelengths, Abs2hv[bin_no,767*2,:]-0.0, "g--", label="2hv %.3fs" % delays[767*2])
    plt.plot(wavelengths, Abs1hv[bin_no,-1,:]-0.0, "r-", label="1hv %.3fs" % delays[-1])
    plt.plot(wavelengths, Abs2hv[bin_no,-1,:]-0.0, "r--", label="2hv %.3fs" % delays[-1])
    plt.ylim(-0.004,0.004)
    plt.xlim(470,750)
    plt.grid(linestyle="dashed")
    plt.legend(loc="upper right")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("\u0394A (550nm)")
    plt.show()



plt.figure(dpi=300)
plt.title(directory.replace("OCP", "run no"))
for bin_no in range(periods.shape[0]):
    plt.plot(delays, Abs1hv[bin_no,:,i_550]-bin_no*0.005, label=str(periods[bin_no])+"s")
plt.ylim(-0.045,0.005)
plt.grid(linestyle="dashed")
plt.legend(loc="upper right")
plt.yticks(np.arange(-0.04, 0.005, 0.005))
plt.xlabel("Time (s)")
plt.ylabel("\u0394A (550nm) 1hv")
plt.show()

plt.figure(dpi=300)
plt.title(directory.replace("OCP", "run no"))
for bin_no in range(periods.shape[0]):
    plt.plot(delays, Abs2hv[bin_no,:,i_550]-bin_no*0.005, label=str(periods[bin_no])+"s")
plt.ylim(-0.045,0.005)
plt.grid(linestyle="dashed")
plt.legend(loc="upper right")
plt.yticks(np.arange(-0.04, 0.005, 0.005))
plt.xlabel("Time (s)")
plt.ylabel("\u0394A (550nm) 2hv")
plt.show()





