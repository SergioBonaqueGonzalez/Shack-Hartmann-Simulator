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

modos=36;

c(1:15)=rand(15,1);

c(16:modos)=0.1*rand(length(16:modos),1);

[ZerRecuperados]=SimAb(modos,c,1025,1.471e-6,0.780,41,1e-3,0,1e-8,16,0);





