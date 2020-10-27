# Generalized Power Law Detector
Software package for the detection and association of marine mammal calls in acoustic data. Project still in process, expect changes in coming weeks.

## Usage
Documentation describing the purpose, inputs, and outputs are included in the comments at the beginning of each Matlab function. The script input_vals.m allows for adjustment of input parameters and storage of outputs for each execution of GPL.m. 
The .wav files included in the repository were recorded at 16,384 Hz on two arrays of five hydrophones each in the Chukchi Sea in the Arctic (Berchok et al., 2016). The calls contained therein originate from bearded seals or bowhead whales. 

## Input Parameters
Reccomended input parameters according to Helble et al. (2012) are gamma = 1; v1 = 1; v2 = 2; eta_thresh = 2.62 * 10^-4; eta_noise = 2.07 * 10^-5; t_min = 0.35;

## Status
Detailed development records can be read in dev_notes.txt. As of 14 October 2020, all functions in the repository are working as expected, and an associator program is the next major development step. 

## Built With
Matlab, including the Signal Processing Toolbox

## Example Figures
![Figure 1](https://i.imgur.com/ngf3Apl.png)
Figure 1: Upper plot shows spectrogram before manipulation with noise and detected signal bounds overlaid. Lower plot shows spectrogram after main spectral content unit of the signal masked onto zero spectrogram with time and frequency bounds overlaid in a different style. Notice that the suggested threshold values did not lead to the detection of the long, faint signal visible here, indicating that the choice of values may be geared more towards minimizing false positive than minimizing false negative detections than is desired in certain applications. Notice also that the vertical streaks from 80-85 sec are not detected as signals due to the minimum signal duration threshold parameter t_min.

## Software Flowchart
![Imgur Image](https://imgur.com/57ZuF7t)
Figure 2: Flowchart depicting all programs and primary data structures
