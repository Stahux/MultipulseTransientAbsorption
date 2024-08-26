# What it is?
This is working code used for multipulse transient absorption experiments, performed in 2024 at Max-Planck Institute for Medical Research, Heidelberg.

# Concept
There are 4 essential scripts used in these experiments.

1) arduino_timing_code.ino is C code running on ATMEL MCU, and generates temporally well defined TTL trigger pulse sequence used to generate laser diode pulses used to excite the sample, to trigger camera acquisition and to switch on and off LED that serves as a probing light source. Typical microcontroller like this (commonly found in Arduino boards) have clock frequency of 16MHz and quartz oscillator. It yields 62.5ns clock period, that will determine timing precision. This timing precision is totally suficient for μs and ms -resolved experiments, and more precise laboratory timing generators are overkill for this task. This script uses MCU hardware counters to generate on-demand TTL pulse seqences desired for these experiments.

2) run_sequence.py is Python script running on regular PC. It manages experimental sequence on higher level, that means moving XY stage (so the sample is probed in the fresh spot between experimental repetitions), initiating TTL pulse sequence (by calling MCU program through RS232 and providing it with pulse sequence parameters) and logging experiment progress. It does not read data from camera (although it may be implemented in the future, however in my case it was just faster and easier to use Andor Solis program running in parallel).

3) preprocess.py is Python script intended to check and reduce obtained data after the experiment is completed (or can be also done using partial data). User needs to specify directory with *.sif files written by Andor Solis software. Each file represents one experimental cycle - or repetition in other words, captured after a series of triggers. Each of these files contain probe intensity defined on temporal (delay) and spectral (wavelength) grid, and some information about camera/spectrometer status. For more details please refer to Andor Solis documentation. *.sif files used as input to this script should be numbered. Script will read all of them into memory, sort them, project data into more coarse wavelength grid, check if sequence is correct based on laser pulse scattering spikes imprinted in the data (and warn user if there is something unexpected), caulculate probe intensity median (using wavelength at which no sample-related signal is present - determined based on previous experiments), filter the data and reject outlier kinetics (that occur because of the small air bubbles in the sample accumulating over time), plot statistics for diagnostics purposes, calculate transient absorbance and save reduced data.

4) plot_data.py is Python script intended to be run after preprocess.py succesfully completed data reduction and processing (that is a slow process). It performs some final steps of data interpretation (like setting proper transient absoorbance signal value before the pulse - so that it is always set to "zero" before the first laser pulse - note that probing does not always cover both laser pulses due to probe light effect on the sample kinetics). Then it plots data in a "pretty way" intended to be used for the presentation/publication.

Note that these scripts are not designed to run out-of-box, they have been written for the specific hardware that I utilized. However, they can be adapted and used to develop similar setups, and I will be happy to help you if you want do develop you own variation of this setup using this code.

# Timing diagram
Experimental sequence implemented by run_sequence.py is depicited below:
![image](https://github.com/user-attachments/assets/41430780-ceae-4ce5-9694-08729f17b02b)
One trigger pulse sequence implemented by arduino_timing_code.ino is depicited below:
![image](https://github.com/user-attachments/assets/78279f44-2783-402f-8bbf-fff3d8bf7634)

# Final comments
Note that these scripts are not designed to run out-of-box, they have been written for the specific hardware that I utilized. However, they can be adapted and used to develop similar setups, and I will be happy to help you if you want do develop you own variation of this setup using this code.
