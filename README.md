# Shack-Hartmann-Simulator in MATLAB

This program was created by Sergio Bonaque-Gonzalez
sergio.bonaque@um.es

It is a simple shack-Hartmann simulator which include the following features;
Customization of:
    Resolution of CCD
    PixelSize
    Wavelength
    Number of microlenses
    Focal of microlenses
Include the option of consider propagation between microlenses array and CCD. The script automatically calcules if Fresnel or Fraunhofer aproximation is necessary
Cuantization 

To Do:
- I like to program in a very easy to understand way. However it is not computationally optimal.
- The script is very commented but in Spanish. If someone is very interested and wants to improve it, a further translation to english is possible.
- Noise it is not considered, but it will be included in the near future.
- The centroid estimation algorithm is basically the one from MATLAB. If noise it is include, the algorithm of centroid estimation should also be improved.
- For testing purposes, in this moment the script only allows to use configurations of pixel and number of microlenses which allow an simmetric distribution of microlesenses and pupil. Now it can be removed, but it is something to do.
- Include separation between microlenses

Syntax:
[ZerRecuperados]=SimAb(modos,c,resolucion,TamPixel,LAMBDA,NLentes,focalML,Propagacion,factor,bits,pintar)

modos= zernike modes (ej. modos=10)
c = vector with the zernike coefficients in Noll notation (ej. c=rand(10,1)
resolucion= resolution of CCD (only square and odd CCDs are considered) (ej. 1024)
TamPixel=Tamaño del pixel en la CCD (ej. 1.471e-6)
LAMBDA= Wavelength in microns (ej. 0.780)
NLentes= Number of microlenses in a row (it is supposed square microlenses array (ej.41))
focalML= Microlenses focal in meters.
Propagacion= Flag that indicates if propagation between microlenses and CCD should be taken into account: =0 no propagation =1 propagation
factor= The value of zernike coeficients used for characterization of the sensor. ej=1e-8
bits= bits used for quantization of CCD
pintar= Flag which indicates if figures should be painted or not 1=yes, 2=no

Example of use:
modos=36;
c(1:15)=rand(15,1);
c(16:modos)=0.1*rand(length(16:modos),1);
[ZerRecuperados]=SimAb(modos,c,2025,1.471e-6,0.780,41,1e-3,0,1e-8,16,0);
