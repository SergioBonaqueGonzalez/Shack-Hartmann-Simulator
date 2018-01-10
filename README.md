# Shack-Hartmann-Simulator in MATLAB


This program was created by Sergio Bonaque-Gonzalez, Optical Engineer & Juan Trujillo-Sevilla, Electronic Engineer

sbonaque@wooptix.com    trujillo@wooptix.com

www.linkedin.com/in/sergiobonaque


We are open to include any request or any contribution.


- The input is a phase represented as a vector of zernike coefficients in the Noll notation 
- The output is a vector of recovered Zernike coefficients through a realistic Shack-Hartmann sensor.

This simulator is in an advanced state of development. It can be useful in Shack-Hartmann design.
It includes the followin items:

    Resolution of CCD. 
    
    PixelSize.
    
    Wavelength.
    
    Number of microlenses.
    
    Focal of microlenses.
    
    Include the option of consider propagation between microlenses array and CCD. 
    
    Cuantization of CCD signal.
    
    Different microlenses geometry
    
    Different methods for centroid calculation
    
    Microlenses can have their own aberrations.
    
    Centroid calculations can be made taking into account crosstalk or not.

    Field distortion can be applied to each microlens
    
    Vignetting of microlenses can be considered
    
    Magnitude of the object (i.e. a star) can be considered.
    
    Pupil diameter
    
    Exposure time
    
    Bandwidth
    
    Quantum efficiency of CCD
    
    PhotonNoise
    
    Well Capacity of CCD
    
    Read Noise of CCD
    
    Dark Current of CCD
    
    Read out noise of CCD
    
    It calculates a variety of quality metrics as PSF, MTF, Sthrel ratio, Convolution with an extended object and analisys of residual error.
    
First, it will calculate and store the recuperation matrix for a certain configuration. This step could take a few minutes. If the configuration is changed, this must be done again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SYNTAX: The simulator is called as a function of the form:

* [RecoveredZern,SH,ML]=SHsimulator(modes,ZerValues,resolution,PixelSize,LAMBDA,Lenses,focal,Prop,factor,bits,...
    MLGeometry,MLCentroidMethod,MLSharedArea,MLfieldDistortion,MLvignetting,ObjectMagnitude,PupilDiameter,...
    Exposuretime,BandWidth,QE,PhotonNoise,WellCapacity,ReadNoise,DarkCurrent,ReadOutNoise,paint)

where:

    modes= number of Zernike modes inteded to be recovered (ej. modes=10)

    ZerValues = vector that contains the Zernike coefficients of the incoming phase in Noll notation. (ej. ZerValues=rand(10,1))

	resolution=resolution of CCD (Only square and odd CCDs are suported) (ej. resolution=1024)

	PixelSize= Pixel size of the CCD (ej. PixelSize=1.471e-6)

	LAMBDA= wavelength of the simulations in microns (ej. LAMBDA=0.780)

	Lenses= number of microlenses in a row in the microlenses array (ej. Lenses=41)

	ML.focal= Focal length of the microlenses in meters. (ej.ML.focal=1e-3)

	Prop= flag that indicates if propagation should be considered between microlenses and CCD =0 NO propagation; =1 Propagation. WARNING: if no propagation exist, microlenses are considered as pure binary amplitude objects.

	factor= the value with which zernize coefficients are characterized through the system. (ej factor=1e-8)

	bits= number of bits of the CCD (ej bits=16)

    MLGeometry= Geometry of the microlenses array.
        1=  perfect spherical lenses in square configuration. Amplitude
            outside microlenses pupil is zero, simulating a lenslet with opacities between microlenses. Defocus term of microlenses is
            the one that maximizing Strhel Ratio for the user-defined focal
            length
        2=  perfect spherical lenses in square configuration. Amplitude
            outside microlenses pupil is one, but phase between microlenses
            is plane(=0).It simulates a lenslet where there exist no
            opacities and areas between microlenses are transpartent but
            without phase.Defocus term of microlenses is
            the one that maximizing Strhel Ratio for the user-defined focal
            length.
        3=  perfect square lenses without space between microlenses.Defocus term of microlenses is
            the one that maximizing Strhel Ratio for the user-defined focal
            length.
    
    MLCentroidMethod= method used for calculating the centroid. (see
    'CentroidCalculation.mat' function for more information
        1= classical way.
        2= displacements from the maximum.
        3= normalization of the peak
        4= setting a threshold in the average of the border, plus 3 times the standard deviation. After that, the centroid is calculated.
        5= setting a threshold in the average of the whole image. Useful
            when borders have no information.
    
    MLAberration= flag that indicates if microlenses has any aberration as,
    for example, spherical aberration or anyone. (0= No, 1= Yes (more realistic)).
    The exact aberrations value can be set in the 'createSH.mat'. As
    default, 0.1 microns of Spherical aberration is considered.

    MLSharedArea= flag that indicates if area behind each microlens take
    into account energy from surrounding microlenses in order to calculate
    centroid (0= No, 1= Yes (more realistic)).

    MLfieldDistortion= flag. if ==1 a function applies Seidel field distortion to the image produced by each microlens in the lenslet.
    (0= No, 1= Yes (more realistic)).

    MLvignetting= flag that indicates if vignetting of microlenses should
    be incorporated to simulations. (0= No, 1= Yes (more realistic)).
	
    ObjectMagnitude= Magnitude of the observed object (i. e. a star)
    (i.e. ObjectMagnitude=0)
    
    PupilDiameter= pupil diameter of the optical system in meters. In a telescope it
    is the diameter of the telescope (i.e. PupilDiameter=4.2).
    
    Exposuretime= exposure time of the system in seconds (i.e.
    Exposuretime=1e-2)
    
    BandWidth= bandwith of the optical system. (i.e. 150). Tipical bandwith
    of filters used in astronomy:
                        filter='U'        BandWidth=54;
                        filter='B'        BandWidth=97;
                        filter='V'        BandWidth=88;
                        filter='R'        BandWidth=147;
                        filter='I'        BandWidth=150;
                        filter='J'        BandWidth=202;
                        filter='H'        BandWidth=368;
                        filter='K'        BandWidth=511;
                        filter='g'        BandWidth=73;
                        filter='r'        BandWidth=94;
                        filter='i'        BandWidth=126;
                        filter='z'        BandWidth=118;
    
    QE= quantum efficiency of the CCD at the used wavelength. (i.e. QE=0.90)
    
    PhotonNoise= Flag. if =1 photon noise will be taken into account. (0=
    NO, 1= Yes (more realistic)).
    
    WellCapacity= Photons well capacity of CCD. It should be scpecified by
    the manufacturer. (i.e. WellCapacity=18000).
    
    ReadNoise=flag that indicates if ReadNoise should be incorporated. (0= No, 1= Yes (more realistic)).
    
    DarkCurrent= Dark current of the CCD in e-/pixel/s. It should be scpecified by
    the manufacturer. (i.e. DarkCurrent=0.01).
    
    ReadOutNoise= Read Noise of the detector in RMS. It should be scpecified by
    the manufacturer. (i.e. ReadOutNoise=8).
        
    paint= dummy flag which indicates if all figures should be painted 1=YES, 2=NO
    


If you run the script directly, you will obtain an example with default values, with the followin outputs:


SHsimulator
Function called without values. Using values by defect. Read documentation to include your own values.
Zernike matrix already exists.
Calculations performed for 16 bits. 
Zernike calibration already exists for this configuration.
Strhel ratio of residual = 0.12691 (ideal =>0.8)
Difference between Strhel ratio of the incoming phase and the recovered one is -0.00111
Difference between the RMS of the incoming phase (1.98537) and the recovered one (1.75808) is = 0.23083 (11.45 of error in percentage).


![My image1](/imgs/figure1.jpg)   
![My image2](/imgs/figure2.jpg)  
![My image3](/imgs/figure3.jpg)  
![My image4](/imgs/figure4.jpg)  
![My image5](/imgs/figure5.jpg)  
![My image6](/imgs/figure6.jpg)  
![My image7](/imgs/figure7.jpg)  
![My image8](/imgs/figure8.jpg)  
![My image9](/imgs/figure9.jpg)  





