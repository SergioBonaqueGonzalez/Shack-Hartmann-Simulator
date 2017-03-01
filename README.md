# Shack-Hartmann-Simulator in MATLAB


This program was created by Sergio Bonaque-Gonzalez, Optical Engineer.

sergio.bonaque@um.es

www.linkedin.com/in/sergiobonaque


I am open to include any request or any contribution.


- The input is a phase represented as a vector of zernike coefficients in the Noll notation 
- The output is a vector of recovered Zernike coefficients through an almost ideal Shack-Hartmann sensor.

In its actual state, it can be used, for example, to test centroid algorithms or new methodology. I hope in the near future will be a very complete simulator.
At the moment it is a shack-Hartmann simulator which includes the following features;
Customization of:

    Resolution of CCD. 
    
    PixelSize.
    
    Wavelength.
    
    Number of microlenses.
    
    Focal of microlenses.
    
    Include the option of consider propagation between microlenses array and CCD. 
    
    The script automatically calcules if Fresnel or Fraunhofer aproximation is necessary.
    
    Cuantization of CCD signal.
    
    It calculates a variety of quality metrics as PSF, MTF, Sthrel ratio, Convolution with an extended object and analisys of residual error.
    
First, it will calculate and store the recuperation matrix for a certain configuration. This step could take a few minutes. If the configuration is changed, this must be done again.




To Do:
- I like to program in a very easy-to-understand way. However it is not computationally optimal.
- Noise it is not considered, but it will be included in the near future.
- The centroid estimation algorithm is basically the one from MATLAB. If noise it is include, the algorithm of centroid estimation should also be improved.
- For testing purposes, in this moment the script only allows to use configurations of pixel and number of microlenses which allow a simmetric distribution of microlesenses and pupil (odd resolution). Now it can be removed, but it is something to do.
- Include separation between microlenses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SYNTAX: The simulator is called as a function of the form:

* [ZerRecuperados]=SimAb(modos,c,resolucion,TamPixel,LAMBDA,NLentes,focalML,Propagacion,factor,bits,pintar)

where:

    modos= how many zernike modes will be considered (ej. modos=10)

    c = vector with the zernike coefficients in Noll notation (ej. c=rand(10,1)

    resolucion= resolution of CCD (only square and odd CCDs are considered) (ej. 1024)

    TamPixel= Pixe size of the CCD in meters (ej. 1.471e-6)

    LAMBDA= Wavelength in microns (ej. 0.780)

    NLentes= Number of microlenses in a row (it is supposed square microlenses array (ej.41))

    focalML= focal length of microlenses  in meters.

    Propagacion= Flag that indicates if propagation between microlenses and CCD should be taken into account: 0=NO ; 1=YES

    factor= The value of zernike coeficients in meters used for characterization of the sensor. ej=1e-8

    bits= bits used for quantization of CCD

    pintar= Flag which indicates if all figures should be painted or not 0=NO ; 1= YES


Example of use:
It will run the simulator for a phase defined by 36 zernike modes. The CCD has 413x413 pixels with a pixel size of 1.471e-6. The wavelength is 0.780 microns. The microlenses array is 15x15. The focal of each microlent is 1e-3 meters. Propagation between microlenses and CCD will be not considered. Each Zernike will be characterized with a value of 1e-8 meters. The CCD has 16 bits and all outputs will be painted.

modos=36;

c(1:15)=rand(15,1);

c(16:modos)=0.1*rand(length(16:modos),1);

[ZerRecuperados]=SimAb(modos,c,413,1.471e-6,0.780,15,1e-3,0,1e-8,16,1);

Zernike matrix not available. Calculating...

Sampling of microlenses is good enough.

Propagation between microlenses and CCD is not considered

Calculations performed for 16 bits. 

Zernike vector is not available. Calculating new Zernike vector...

Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.

lsqr converged at iteration 57 to a solution with relative residual 0.034.

Strhel ratio of residual = 0.91709

Difference between Strhel ratio of the incoming phase and the recovered one is -0.00019

Difference between the RMS of the incoming phase (1.93317) and the recovered one (1.90105) is = 0.04848 (1.66 of error in percentage).

The size of the pupil is 0.00075 meters

Elapsed time is 29.992629 seconds.








