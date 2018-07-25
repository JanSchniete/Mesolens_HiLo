# Mesolens_HiLo
MATLAB script for end-user HiLo microscopy and mesoscopy. Currently tuned for Mesolens data but easily adjustable. Few homebuilt functions are
supplied with the script. This script is meant to be used with laser speckle illumination HiLo microscopy but can be used as a starting point
to write a script for grid illumination HiLo microscopy.
The three .m files bpgauss.m, lpgauss.m and hpgauss.m are filter creation functions used in the actual scrip hilostack.m
How to use:
1. Upon running the script hilostack.m you will be asked several times to select various files.

    1.1 First you can select a correction file for flatness of field. This is optional, if you press 'cancel', the correction defaults to unity.

    1.2 Second you have to select ALL uniform illumination files you wish to process.

    1.3 Third you have to select ALL speckle illumination files for processing

NOTE: Number of files of 1.2 and 1.3 have to be the same (these are image pairs)

2. Next you have to specify the HiLo parameters used for processing

    2.1 First is the parameter for optical sectioning thickness. This should be done in integer numbers starting from 1. The actual thickness 
    will depend on your image parameters as well as the speckle pattern forming in the sample.

    2.2 Second the scale of your images in pixels/micron.

3. The script will run now until all HiLo related processing is complete and will evaluate the ratio of high to low spatial frequencies at the
    cut-off frequency between low and high pass filter. You can now pick a scaling or choose an individual one, this should
    be close to 1 as the images have been normalised. However, it is best to try different values and evaluate which gives the best results.
    
NOTE: This code makes use of the gpuArray function of MATLAB, if you receive an out-of-memory error for either RAM or VRAM, please contact me.
I have several versions of the code capable of running on weaker machines.
