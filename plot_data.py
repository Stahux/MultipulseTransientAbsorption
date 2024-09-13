# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 17:19:09 2024

@author: Stanisław Niziński
"""

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import copy as cp

#########################################################################
############THIS IS 2nd STEP PLOTTING/FINAL PROCESSING CODE##############
################ it is meant to use reduced data after the ##############
################ preprocess.py script to produce final     ##############
################ plots and final data files                ##############
################ It is not meant                           ##############
################ to be fast or elegant. It is meant to     ##############
################ do its job correctly. It is always work   ##############
################ in progress and tailored for this         ##############
################ specific experiment.                      ##############
#########################################################################

################### LOAD DATA ###################
#################################################
directory = "exemplary_data" #here specify dataset that should be processed
#################################################
buffer_directory = ""#../12-05-2024 buffer 47" #specify buffer-only dataset for subtraction
#################################################
#################################################
#################################################


#below set if use all delays or only those for where good correction is done
nkins = 10 #use 10 for 2hv proteins and 8 for 1hv proteins
spacing = 0.002 #for multi-kinetic graphs - vertical offset

#save drop this code to dir too to know later what was used to process this data
file_code = open("plot_data.py", "r")
file_codeout = open(directory+"\\plot_data_copy.py", "w") #code copy for logging
code_text = file_code.read()
file_codeout.write(code_text)
file_code.close()
file_codeout.close()


#load important data files
delays = np.load(directory+"\\delays.npy")
wavelengths = np.load(directory+"\\wavelengths.npy")
periods = np.load(directory+"\\periods.npy")
periodpattern = np.tile(periods, 90) #like in the sequence code

Abs1hv = np.load(directory+"\\Abs1hv.npy")  
Abs2hv = np.load(directory+"\\Abs2hv.npy")    


#just to see that you use right - corrected grid
print("Longest probe wavelength in grid is: ", wavelengths[-1])

#enable to subtract equivalent buffer dataset - any drifts not related directly to OCP
if(buffer_directory != ""):
    buffer_delays = np.load(buffer_directory+"\\delays.npy")
    buffer_wavelengths = np.load(buffer_directory+"\\wavelengths.npy")
    buffer_periods = np.load(buffer_directory+"\\periods.npy")
    buffer_periodpattern = np.tile(periods, 90)
    
    buffer_Abs1hv = np.load(buffer_directory+"\\Abs1hv.npy")  
    buffer_Abs2hv = np.load(buffer_directory+"\\Abs2hv.npy")        

    if(not(buffer_Abs1hv.shape == Abs1hv.shape)):
        raise Exception("Delays of buffer dataset do not match!")    
    if(not((buffer_delays == delays).all())):
        raise Exception("Delays of buffer dataset do not match!")
    if(not((buffer_wavelengths == wavelengths).all())):
        raise Exception("Wavelengths of buffer dataset do not match!")    
    if(not((buffer_periods == periods).all())):
        raise Exception("Periods of buffer dataset do not match!")    
        
    Abs1hv = Abs1hv - buffer_Abs1hv
    Abs2hv = Abs2hv - buffer_Abs2hv


#average 1 pulse kinetics to get one of high quality:
Abs1hv_avg = np.average(Abs1hv,axis=0)

Abs2hv_corr = Abs2hv.copy()

dshape = Abs1hv[0].shape
delays = delays-2.490 # set proper time zero, 10ms is probably camera lag minus laser lag


#below code sets proper zero absorbance before the pulse (uses 1hv kinetic if possible)
alignment_time = -0.1 #zeroing/aligning at -100ms
alignment_i = np.argmin(np.abs(delays-alignment_time))

if(True): #set zero for ref kinetics
    Abs1hv_zero = np.average(Abs1hv_avg[alignment_i-5:alignment_i+6,:], axis=0)
    Abs1hv_avg[:,:] -= np.tile(Abs1hv_zero,(dshape[0],1))    
    for bin_no in range(periods.shape[0]):
        Abs1hv_zero = np.average(Abs1hv[bin_no,alignment_i-5:alignment_i+6,:], axis=0)
        Abs1hv[bin_no,:,:] -= np.tile(Abs1hv_zero,(dshape[0],1))

if(True): #set zero according to ref kinetic...
    for bin_no in range(periods.shape[0]):
        if(periods[bin_no] < 1.5): #no gluing needed there, because both pulses are in the window
            alignment_delay = alignment_time-periods[bin_no]
            alignment_i_tmp = np.argmin(np.abs(delays-alignment_delay))
            zero_value = np.average(Abs2hv_corr[bin_no,alignment_i_tmp-5:alignment_i_tmp+6,:], axis=0)
            Abs2hv_corr[bin_no] -= zero_value
        elif(periods[bin_no] < 15): #period must be lower than single pulse kinetic duration, otherwise no gluing posssible
            #firstly find signal in 1 pulse kinetic at delay of 2nd pulse - alignment_time
            ref_i = np.argmin(np.abs(delays-(periods[bin_no]-alignment_time)))
            ref_value = np.average(Abs1hv_avg[ref_i-5:ref_i+6,:], axis=0)
            #and respective time point in 2hv kinetic
            zero_value = np.average(Abs2hv_corr[bin_no,alignment_i-5:alignment_i+6,:], axis=0)
            #and now ensure that zero_value is set to ref_value in new kinetic
            Abs2hv_corr[bin_no] -= np.tile(zero_value - ref_value, (dshape[0],1))
        else: #unfortunately cannot correct above 15s properly. 
            #hoping that at least at some wavelength things relaxed before 2nd pulse
            #if not, then just dont use these data parts!
            zero_value = np.average(Abs2hv_corr[bin_no,alignment_i-5:alignment_i+6,:], axis=0)
            Abs2hv_corr[bin_no] -= np.tile(zero_value, (dshape[0],1))      
            
if(True): #subtract zero in a regular way - produces hard-zeroed array.
    for bin_no in range(periods.shape[0]):
        Abs2hv_zero = np.average(Abs2hv[bin_no,5:15,:], axis=0)
        Abs2hv[bin_no,:,:] -= np.tile(Abs2hv_zero,(dshape[0],1))



if(False): #bucket avg filter - enable to use smooth average in the temporal regime
    boff = 1 #how many left and right to take into bucket, 2xboff+1 will go to one bucket
    
    points_no = Abs1hv.shape[1]
    
    Abs1hv_new = cp.deepcopy(Abs1hv)*np.nan
    for i in range(0+boff,points_no-boff):
        Abs1hv_new[:,i,:] = np.average(Abs1hv[:,i-boff:1+i+boff,:],axis=1)
    Abs1hv = Abs1hv_new[:,0+boff:points_no-boff,:]

    Abs1hv_new = cp.deepcopy(Abs1hv_avg)*np.nan
    for i in range(0+boff,points_no-boff):
        Abs1hv_new[i,:] = np.average(Abs1hv_avg[i-boff:1+i+boff,:],axis=0)
    Abs1hv_avg = Abs1hv_new[0+boff:points_no-boff,:]

    Abs2hv_new = cp.deepcopy(Abs2hv)*np.nan
    for i in range(0+boff,points_no-boff):
        Abs2hv_new[:,i,:] = np.average(Abs2hv[:,i-boff:1+i+boff,:],axis=1)
    Abs2hv = Abs2hv_new[:,0+boff:points_no-boff,:]

    Abs2hv_new = cp.deepcopy(Abs2hv_corr)*np.nan
    for i in range(0+boff,points_no-boff):
        Abs2hv_new[:,i,:] = np.average(Abs2hv_corr[:,i-boff:1+i+boff,:],axis=1)
    Abs2hv_corr = Abs2hv_new[:,0+boff:points_no-boff,:]

    delays = delays[0+boff:points_no-boff]
    


#cutting
t_cut = np.argmin(np.abs(delays+1.5)) #cut at -1s

delays = delays[t_cut:]
Abs1hv = Abs1hv[:,t_cut:,:]
Abs2hv = Abs2hv[:,t_cut:,:]
Abs1hv_avg = Abs1hv_avg[t_cut:,:]
Abs2hv_corr = Abs2hv_corr[:,t_cut:,:]
dshape = Abs1hv[0].shape


####################!!!!!!!!!!!!!!!!!!!!!!!!!!!!#######################
#keep line below uncommented to use corrected kinetic instead of hard-zeroed!
Abs2hv = Abs2hv_corr



#################################################
################### PLOT DATA ###################
#################################################

#get coords funcs
i_d = lambda delay: np.argmin(np.abs(delays-delay))
i_w = lambda wavelength: np.argmin(np.abs(wavelengths-wavelength))


#calculate product vs delay curve
dA550list = []
for bin_no in range(nkins):
    #delay = periods[bin_no]
    dA550 = np.average(Abs2hv[bin_no,i_d(5)-10:i_d(5)+11,i_w(585)-2:i_w(585)+3])
    dA550list.append(dA550)
dA550list = np.array(dA550list)

plt.figure(dpi=300)
#plt.title(directory)
plt.plot(periods[:nkins], dA550list, "o--", label="")
#plt.ylim(-0.045,0.005)
plt.grid(linestyle="dashed")
#plt.legend(loc="upper left")
plt.xlabel("Delay between pulses (s)")
plt.xscale("log")
plt.ylabel("Final product - \u0394A(580nm,5s)")
plt.savefig(directory+"\\product_style2hv.png")
plt.show()
np.save(directory+"\\product_curve", dA550list)


for bin_no in range(periods.shape[0]):
    
    plt.figure(dpi=300)
    #plt.title(directory.replace("OCP", "run no")+", delay "+str(periods[bin_no])+"s")
    plt.plot(delays, np.average(Abs1hv_avg[:,i_w(550)-1:i_w(550)+2],axis=1)-0.0, label=r"$1h\nu$")
    plt.plot(delays, np.average(Abs2hv[bin_no,:,i_w(550)-1:i_w(550)+2],axis=1)-0.0, label=r"$2h\nu$")
    #plt.ylim(-0.01,0.01)
    plt.xlim(-1.5,17.5)
    plt.grid(linestyle="dashed")
    plt.legend(loc="lower right")
    plt.text(0.98, 0.95, r'$t_{delay}='+str(periods[bin_no])+'s$', horizontalalignment='right',
         verticalalignment='top', transform=plt.gca().transAxes, fontsize=14)
    plt.xlabel("Time (s)")
    plt.ylabel("\u0394A (550nm)")
    plt.savefig(directory+"\\kinetics550_"+str(periods[bin_no])+"s"+"_style2hv.png")
    plt.show()
    
    plt.figure(dpi=300)
    #plt.title(directory.replace("OCP", "run no")+", delay "+str(periods[bin_no])+"s")
    plt.plot(delays, np.average(Abs1hv_avg[:,i_w(585)-1:i_w(585)+2],axis=1)-0.0, label=r"$1h\nu$")
    plt.plot(delays, np.average(Abs2hv[bin_no,:,i_w(585)-1:i_w(585)+2],axis=1)-0.0, label=r"$2h\nu$")
    #plt.ylim(-0.01,0.01)
    plt.xlim(-1.5,17.5)
    plt.grid(linestyle="dashed")
    plt.legend(loc="lower right")
    plt.text(0.98, 0.95, r'$t_{delay}='+str(periods[bin_no])+'s$', horizontalalignment='right',
         verticalalignment='top', transform=plt.gca().transAxes, fontsize=14)
    plt.xlabel("Time (s)")
    plt.ylabel("\u0394A (585nm)")
    plt.savefig(directory+"\\kinetics585_"+str(periods[bin_no])+"s"+"_style2hv.png")
    plt.show()    
    
    plt.figure(dpi=300)
    #plt.title(directory.replace("OCP", "run no")+", delay "+str(periods[bin_no])+"s")
    #plt.plot(wavelengths, Abs1hv_avg[579,:]-0.0, "b--", label="1hv %.3fs" % delays[579])
    #plt.plot(wavelengths, Abs2hv[bin_no,579,:]-0.0, "b-", label="2hv %.3fs" % delays[579])
    #plt.plot(wavelengths, Abs1hv_avg[581,:]-0.0, "b--", label="1hv %.3fs" % delays[581])
    #plt.plot(wavelengths, Abs2hv[bin_no,581,:]-0.0, "b-", label="2hv %.3fs" % delays[581])
    plt.plot(wavelengths, Abs1hv_avg[i_d(0.5),:]-0.0, "b--", label=r"$1h\nu,t=%.2fs$" % 0.5)
    plt.plot(wavelengths, Abs2hv[bin_no,i_d(0.5),:]-0.0, "b-", label=r"$2h\nu,t=%.2fs$" % 0.5)
    #plt.plot(wavelengths, Abs1hv_avg[767,:]-0.0, "c--", label="1hv %.3fs" % delays[767])
    #plt.plot(wavelengths, Abs2hv[bin_no,767,:]-0.0, "c-", label="2hv %.3fs" % delays[767])
    plt.plot(wavelengths, Abs1hv_avg[i_d(5),:]-0.0, "r--", label=r"$1h\nu,t=%.2fs$" % 5)
    plt.plot(wavelengths, Abs2hv[bin_no,i_d(5),:]-0.0, "r-", label=r"$2h\nu,t=%.2fs$" % 5)
    #plt.plot(wavelengths, Abs1hv_avg[2870,:]-0.0, "r--", label="1hv %.3fs" % delays[2870])
    #plt.plot(wavelengths, Abs2hv[bin_no,2870,:]-0.0, "r-", label="2hv %.3fs" % delays[2870])
    #plt.plot(wavelengths, Abs1hv_avg[-1,:]-0.0, "g--", label="1hv %.3fs" % delays[-1])
    #plt.plot(wavelengths, Abs2hv[bin_no,-1,:]-0.0, "g-", label="2hv %.3fs" % delays[-1])
    plt.xlim(470,750)
    plt.grid(linestyle="dashed")
    plt.legend(loc="lower right")
    plt.text(0.98, 0.4, r'$t_{delay}='+str(periods[bin_no])+'s$', horizontalalignment='right',
         verticalalignment='top', transform=plt.gca().transAxes, fontsize=14)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("\u0394A")
    plt.savefig(directory+"\\spectra_"+str(periods[bin_no])+"s"+"_style2hv.png")
    plt.show()
    
    """ ##useful for pulsetrain plots
    plt.figure(dpi=300)
    plt.title(directory.replace("OCP", "run no")+", delay "+str(periods[bin_no])+"s")
    plt.plot(wavelengths, Abs1hv_avg[681,:]-0.0, "r--", label="1hv %.3fs" % delays[681])
    plt.plot(wavelengths, Abs2hv[bin_no,681,:]-0.0, "r-", label="2hv %.3fs" % delays[681])
    plt.plot(wavelengths, Abs1hv_avg[681+384,:]-0.0, "y--", label="1hv %.3fs" % delays[681+384])
    plt.plot(wavelengths, Abs2hv[bin_no,681+384,:]-0.0, "y-", label="2hv %.3fs" % delays[681+384])
    plt.plot(wavelengths, Abs1hv_avg[681+2*384,:]-0.0, "g--", label="1hv %.3fs" % delays[681+2*384])
    plt.plot(wavelengths, Abs2hv[bin_no,681+2*384,:]-0.0, "g-", label="2hv %.3fs" % delays[681+2*384])
    plt.plot(wavelengths, Abs1hv_avg[681+3*384,:]-0.0, "c--", label="1hv %.3fs" % delays[681+3*384])
    plt.plot(wavelengths, Abs2hv[bin_no,681+3*384,:]-0.0, "c-", label="2hv %.3fs" % delays[681+3*384])
    plt.plot(wavelengths, Abs1hv_avg[681+4*384,:]-0.0, "b--", label="1hv %.3fs" % delays[681+4*384])
    plt.plot(wavelengths, Abs2hv[bin_no,681+4*384,:]-0.0, "b-", label="2hv %.3fs" % delays[681+4*384])
    plt.xlim(470,750)
    plt.grid(linestyle="dashed")
    plt.legend(loc="lower right")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("\u0394A")
    plt.savefig(directory+"\\spectra_"+str(periods[bin_no])+"s"+"_style2hv.png")
    plt.show()
    """


#below is obsolete because I average all 1hv kinetics!
#plt.figure(dpi=300)
#plt.title(directory.replace("OCP", "run no"))
#for bin_no in range(periods.shape[0]):
#    plt.plot(delays, Abs1hv[bin_no,:,i_w(580)]-bin_no*spacing, label=str(periods[bin_no])+"s")
#plt.ylim(-10*spacing,spacing)
#plt.xlim(-1,6)
#plt.grid(linestyle="dashed")
#plt.legend(loc="upper right")
#plt.yticks(np.arange(-9*spacing, spacing, spacing))
#plt.xlabel("Time (s)")
#plt.ylabel(r"\u0394A (580nm), $1h\nu$")
#plt.savefig(directory+"\\kinetics_refkinetics_style2hv.png")
#plt.show()



plt.figure(dpi=300, figsize=(4, 4.8))
#plt.title(directory.replace("OCP", "run no"))
for bin_no in range(nkins):
    plt.plot(delays, np.average(Abs2hv[bin_no,:,i_w(585)-1:i_w(585)+2],axis=1)-bin_no*spacing, label=str(periods[bin_no])+"s")
plt.ylim(-nkins*spacing,spacing)
plt.xlim(-1,2)
plt.grid(linestyle="dashed")
plt.legend(loc="upper right", title=r'$t_{delay}:$')
plt.yticks(np.arange(-(nkins-1)*spacing, spacing, spacing))
plt.xlabel("Time (s)")
plt.ylabel("\u0394A (585nm), $2h\\nu$")
plt.savefig(directory+"\\kinetics_kinetics_style2hv.png")
plt.show()


plt.figure(dpi=300)
#plt.title(directory.replace("OCP", "run no"))
for bin_no in range(nkins):
    plt.plot(delays, np.average(Abs2hv[bin_no,:,i_w(585)-1:i_w(585)+2],axis=1)-bin_no*spacing, label=str(periods[bin_no])+"s")
plt.ylim(-nkins*spacing,spacing)
plt.xlim(-1,5)
plt.grid(linestyle="dashed")
plt.legend(loc="upper right", title=r'$t_{delay}:$')
plt.yticks(np.arange(-(nkins-1)*spacing, spacing, spacing))
plt.xlabel("Time (s)")
plt.ylabel("\u0394A (585nm), $2h\\nu$")
plt.savefig(directory+"\\kinetics_kinetics_style2hv.png")
plt.show()



#detailed plot of the reference kinetic
plt.figure(dpi=300)
#plt.title(directory.replace("OCP", "run no")+", reference")
plt.plot(wavelengths, Abs1hv_avg[i_d(0.01),:]-0.0, "b-", label=r"$1h\nu,t=%.2fs$" % 0.01)
plt.plot(wavelengths, Abs1hv_avg[i_d(0.03),:]-0.0, "c-", label=r"$1h\nu,t=%.2fs$" % 0.03)
plt.plot(wavelengths, Abs1hv_avg[i_d(0.1),:]-0.0, "g-", label=r"$1h\nu,t=%.2fs$" % 0.1)
plt.plot(wavelengths, Abs1hv_avg[i_d(0.5),:]-0.0, "y-", label=r"$1h\nu,t=%.2fs$" % 0.5)
plt.plot(wavelengths, Abs1hv_avg[i_d(1),:]-0.0, "r-", label=r"$1h\nu,t=%.2fs$" % 1)
plt.plot(wavelengths, Abs1hv_avg[i_d(5),:]-0.0, "m-", label=r"$1h\nu,t=%.2fs$" % 5)
#plt.plot(wavelengths, Abs1hv_avg[-1,:]-0.0, "k-", label=r"$1h\nu,t=%.2fs$" % delays[-1])
#plt.ylim(-0.0025,0.002)
plt.xlim(470,750)
plt.grid(linestyle="dashed")
plt.legend(loc="lower right")
plt.xlabel("Wavelength (nm)")
plt.ylabel("\u0394A")
plt.savefig(directory+"\\reference_style2hv.png")
plt.show()    



