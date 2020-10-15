# Generalized Power Law Detector
Software package for the detection and association of marine mammal calls in acoustic data

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
![][/tlwp_masked.png]
