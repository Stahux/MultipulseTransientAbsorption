# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 11:39:30 2024

@author: Stanisław Niziński
"""


import serial
import time
import numpy as np
import datetime


###############################################################################
################        THIS IS EXPERIMENT RUNNING CODE              ##########
################ It runs experimental sequence, that means           ##########
################ calling arduino to generate timing sequence for     ##########
################ each cycle, moving XY stage in between and          ##########
################ logging. Camera software (Andor Solis) operates     ##########
################ independently in the loop, capturing kinetic        ##########
################ after every trigger pulse and saving each           ##########
################ kinetic to sif file with next number.               ##########
################ Looping Andor Solis is easy, for example one        ##########
################ can set "Kinetic" mode, with "External Start"       ##########
################ and FVB binning, enable autosave with autoincrement.##########
################ Then save a following script as *.pgm file:         ##########
################                                                     ##########
################    counter = 0                                      ##########
################    while counter < 9999999                          ##########
################    	run()                                        ##########
################    	counter = counter + 1                        ##########
################    wend                                             ##########
################                                                     ##########
################ And run it from Andor Solis. It will loop           ##########
################ the camera acquisitions. One must ensure that       ##########
################ camera captures kinetics before next trigger        ##########
################ arrives!                                            ##########
################                                                     ##########
################ Note that script below is not meant                 ##########
################ to be fast or elegant. It is meant to               ##########
################ do its job correctly. It is always work             ##########
################ in progress and tailored for this                   ##########
################ specific experiment.                                ##########
###############################################################################

#here drop some logging data about each executed seq.
file = open("20-08-2024 OCP 79.log", "w") 

#save code to file too
file_code = open("run_sequence.py", "r")
code_text = file_code.read()
file.write(code_text)
file.write("\n\n\n")
file_code.close()

#ad some comment here!
file.write("Sequence type II! OCP R27L Plk ECN notag OCP80 new aliq. diluted and 10x less probe. Usual long seq! Oxxious 520nm and 64us pulse!")
file.write("\n\n")


#open serial port for motion controller
esp = serial.Serial('COM8', 921600, timeout=1)

# Open the serial port for arduino - to send timing commands
ser = serial.Serial('COM5', timeout=1)

mess ="0 0\n" #LED OFF
ser.write(mess.encode())

#this is very lame way of setting each t_delay 3 times
#it is like this because stage needs 3 positions before it moves to next
#delay value (one position with 0 pulses, with 1 pump pulse, 2 one with 2 pulses, 
#and then it will stay in the same place and record dark
period_grid = np.array([0.00300, 0.00300, 0.00300,
                        0.01000, 0.01000, 0.01000, 
                        0.03000, 0.03000, 0.03000, 
                        0.10000, 0.10000, 0.10000,
                        0.30000, 0.30000, 0.30000,
                        1.00000, 1.00000, 1.00000,
                        3.00000, 3.00000, 3.00000,
                        10.0000, 10.0000, 10.0000,
                        30.0000, 30.0000, 30.0000,
                        100.000, 100.000, 100.000])
#so repeat is 270x times - but this is a lot and one can stop sequence earlier 
#and still have a nice dataset
periodlist = np.tile(period_grid, 78)

#set up stage grid params - total points number for the XY stage must be 3N 
hgrid = np.linspace(0,12,13)
vgrid = np.linspace(2.3,10.4,18)

#just make nice arrays of XY positions to iterate over
hh, vv = np.meshgrid(hgrid, vgrid)
#horiz vert value list
poslist = np.zeros([hgrid.shape[0]*vgrid.shape[0],2])
poslist[:,0] = hh.flatten()
poslist[:,1] = vv.flatten()
#and then repeat it a few times to match periodlist size
poslist = np.concatenate([poslist,poslist,poslist,poslist,poslist,poslist,poslist,poslist,poslist,poslist]) #x10 because 10 delays



file.write(str(poslist))
file.write("\n")

file.write(str(periodlist))
file.write("\n")

if(periodlist.shape[0] != poslist.shape[0]):
    raise Exception("Grids have not equal size!!!")


if(len(poslist)/3 - round(len(poslist)/3) != 0.0):
    raise Exception("Stage grid point number must be 3N!!!")

#I used two arrays with different sizes, because before cycle (4) not stage
#movement is executed. Not very elegant, but works.
count_max = 4*int(len(poslist)/3) #add 4th seq. element without stage movement
i = 0 #this is just stage position counter, not always incremented.
count = 0

#give yourself some time to leave and lock the room and return when finished!
time.sleep(3)


while count < count_max:
    
    #ask motor to go to new position
    mesp = "1PA%.5f\r" % poslist[i][1]
    print(esp.write(mesp.encode()))
    mesp = "2PA%.5f\r" % poslist[i][0]
    print(esp.write(mesp.encode()))
    
    p1_old = np.inf
    p2_old = np.inf
    
    #wait for motors
    while True:
        time.sleep(0.1)
        
        esp.write("1TP\r".encode())
        p1 = float(esp.readline().decode())
        esp.write("2TP\r".encode())
        p2 = float(esp.readline().decode())
        
        if(p1 == p1_old and p2 == p2_old):
            print("Movement finished! Arrived to %.5f, %.5f !" % (p2, p1))
            file.write("Movement finished! Arrived to %.5f, %.5f !\n" % (p2, p1))
            break
        else:
            p1_old = p1
            p2_old = p2
    #motors movement done! you can proceed!
    
    print("Now pulse spacing will be %.6f" % periodlist[i])
    
    period_set = 0.001 #period is 1ms so it will be always like this - I will use multiples of it
    laser_pulse_length = 0.000064 #laser pulse length, 6 digit precision will be taken
    #this laser pulse length determines laser pulse energy! 64us is 50mJ/cm2.

    led_1_no = 1 #number of the counter overflow when led on
    camera_1_no = 101 #camera will get pulse in ov. 101 (so 100ms warm up time for LED)
    laser_1_no = int(2.5/period_set)+camera_1_no #laser will get 1st pulse 2.5 after camera trig
    led_1off_no = int(20.5/period_set)+camera_1_no #led off 20.5s after camera trig
    
    laser_2_no = int(periodlist[i]/period_set)+laser_1_no # calc ov. of the second pulse
    camera_2_no = laser_2_no-int(2.5/period_set) #2.5 before laser trig
    led_2_no = camera_2_no-100 #when led on at 2 pulse
    led_2off_no = int(20.5/period_set)+camera_2_no #led off 20.5s after camera trig
    
    
    #so there must be led_xoff_no+1 total_overflows!
    
    total_overflows_1 = led_1off_no+1 #for short seq.
    total_overflows_2 = led_2off_no+1 #for 2 pulse long seq

    #in this seq. idea is to have probe 2.5s before the laser pulse
    #and if one wants it off 7.5s after the pulse, 3832 spectra in andor are needed
    #and if one wants it off 17.5s after the pulse, 7664 spectra in andor are needed
    #and led will be switched on always 10 overflows before camera
    
    #waiting times for the sequence loop in this script
    waiting_time_1 = 10+20 #for ref 1hv and dark
    waiting_time_2 = waiting_time_1 + periodlist[i] + 5 #for 2hv
 
    
    #timestamp of sequence request!    
    print(datetime.datetime.now())
    file.write(str(datetime.datetime.now())+"\n")    
    
    ######CYCLE (1)
    if(count % 4 == 0): #firstly do ref seq. with no laser
        file.write("Seq. type 0!\n")
        
        #setup sequence params for led
        message ="4 %i %i\n" % (led_1_no,led_1off_no)
        ser.write(message.encode())
        print(message)
        file.write(message)          
        
        #setup sequence params for laser 
        message ="1 1 %i\n" % (total_overflows_1+1) #will never happen
        ser.write(message.encode())
        print(message)
        file.write(message)  
    
        #setup sequence params for camera 
        message ="2 1 %i\n" % (camera_1_no)
        ser.write(message.encode())
        print(message)
        file.write(message) 
        
        #start seq command
        message ="3 %.6f %.6f 0.0001 %i\n" % (period_set, laser_pulse_length, total_overflows_1)
        ser.write(message.encode())
        print(message)
        file.write(message)
        
        waiting_time = waiting_time_1
        
    ######CYCLE (2)
    elif(count % 4 == 1): #secondly do regular seq. with 1 laser pulse
        file.write("Seq. type 1!\n")

        #setup sequence params for led
        message ="4 %i %i\n" % (led_1_no,led_1off_no)
        ser.write(message.encode())
        print(message)
        file.write(message)          
        
        #setup sequence params for laser 
        message ="1 1 %i\n" % laser_1_no #only one laser
        ser.write(message.encode())
        print(message)
        file.write(message)  
    
        #setup sequence params for camera 
        message ="2 1 %i\n" % (camera_1_no)
        ser.write(message.encode())
        print(message)
        file.write(message) 
        
        #start seq command
        message ="3 %.6f %.6f 0.0001 %i\n" % (period_set, laser_pulse_length, total_overflows_1)
        ser.write(message.encode())
        print(message)
        file.write(message)        
        
        waiting_time = waiting_time_1
        
    ######CYCLE (3)
    elif(count % 4 == 2): #thirdly do seq. with 2 pulses
        file.write("Seq. type 2!\n")       

        #setup sequence params for led
        message ="4 %i %i\n" % (led_2_no,led_2off_no)
        ser.write(message.encode())
        print(message)
        file.write(message)          
        
        #setup sequence params for laser 
        message ="1 2 %i %i\n" % (laser_1_no,laser_2_no) #two laser pulses!
        ser.write(message.encode())
        print(message)
        file.write(message)  
    
        #setup sequence params for camera 
        message ="2 1 %i\n" % (camera_2_no)
        ser.write(message.encode())
        print(message)
        file.write(message) 
        
        #start seq command
        message ="3 %.6f %.6f 0.0001 %i\n" % (period_set, laser_pulse_length, total_overflows_2)
        ser.write(message.encode())
        print(message)
        file.write(message)     
        
        waiting_time = waiting_time_2
        
    ######CYCLE (4)
    else: #just dark kinetic with no laser no probe
        file.write("Seq. type 3!\n")
        
        #setup sequence params for led
        message ="4 -1 -1\n" #no led there
        ser.write(message.encode())
        print(message)
        file.write(message)          
        
        #setup sequence params for laser 
        message ="1 1 %i\n" % (total_overflows_1+1) #will never happen
        ser.write(message.encode())
        print(message)
        file.write(message)  
    
        #setup sequence params for camera 
        message ="2 1 %i\n" % (camera_1_no)
        ser.write(message.encode())
        print(message)
        file.write(message) 
        
        #start seq command
        message ="3 %.6f %.6f 0.0001 %i\n" % (period_set, laser_pulse_length, total_overflows_1)
        ser.write(message.encode())
        print(message)
        file.write(message)
    
        waiting_time = waiting_time_1
    

    
    # Read and print the serial output - wait long enough for all things to complete
    start_time = time.time()
    while (time.time() - start_time) < waiting_time:  # Wait
        if ser.in_waiting > 0:
            readed = ser.readline()
            print(readed)
            file.write(readed.decode())
        else:
            time.sleep(0.1)
    
    
    print("Finished cycle with count %i and i %i." % (count, i)) 
    file.write("Finished cycle with count %i and i %i.\n" % (count, i))
    
    if(count % 4 != 2): #do not go to other stage point after 3rd seq. element!
        i += 1 #so 3rd and 4th elements of the seq. will be at the same stage position!
    count += 1 #increment loop count

print("Finished measurement with no errors!")
file.write("Finished measurement with no errors!\n")

# Close the connection
ser.close()
esp.close()
#upled.close()
file.close()


###################################################################################
# This script can be killed during the experiment without harm,                   #
# the only harm will be less data. But remember to call then manually file.close()#
#So that log is completed gracefully. Processing script will just load the data   #
#that was saved before killing the script and it will check sequence correctness  #
###################################################################################
