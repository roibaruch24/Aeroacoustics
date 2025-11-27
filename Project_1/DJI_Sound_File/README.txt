%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% README FILE FOR HW1 - AIRCRAFT AEROACOUSTICS - 0860395
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hadar Ben-Gida

The file 'Run35_DJI_white_baseline_B2_CW-calibrated.mat' includes acoustic data measured in the anechoic chamber at the Technion.

Please load the file into MATLAB using:
load('.../Run35_DJI_white_baseline_B2_CW-calibrated.mat')

The .mat file includes several variables:
1. B           - Number of blades (=2)
2. D           - Rotor diameter in [m] (=0.2413m)
3. Fqs         - Sampling rate of the microphone data, in [Hz]
4. MicData     - An array of size (2428166x15) that consists of the recorded microphone data. 	
                 Rows are time samples, and columns correspond to different microphones, starting from theta=0deg and moving to theta=105deg.
5. num_of_mics - Number of microphones exist (=15)
6. RPM.        - Rotational speed of the rotor in rev/min
7. theta       - Azimuth angles of the different microphones. 
                 Each column here is associated with the corresponding column in MicData.



To play the microphone recording and hear the noise of the rotor, you can follow the steps, but first MAKE SURE YOUR SPEAKERS ARE TURNED TO MINIMUM! After all, you don't want to be annoyed by the rotor noise :P

To play the sound for the Microphone above the rotor (theta=0deg), you can enter the following command in MATLAB:
sound(MicData(:,1), Fqs);

To stop the sound from playing, run the following command:
clear sound



