function [RecoveredZern,SH,ML]=SHsimulator(modes,ZerValues,resolution,PixelSize,LAMBDA,Lenses,focal,Prop,factor,bits,...
    MLGeometry,MLCentroidMethod,MLSharedArea,MLfieldDistortion,MLvignetting,ObjectMagnitude,PupilDiameter,...
    Exposuretime,BandWidth,QE,PhotonNoise,WellCapacity,ReadNoise,DarkCurrent,ReadOutNoise,paint)
%{
Created by:
 Sergio Bonaque-Gonzalez. Optical Engineer.
 sergiob@wooptix.com
    &
 Juan Trujillo-Sevilla. Electronic Engineer
  trujillo@wooptix.com
 
    
This program simulates a Shack-Hartmann sensor.

INPUTS
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


IMPORTANT ISSUES AND TO DO:
    - You can incorporate any aberration to microlenses, see documentation
    of "*Definition of the microlenses array" section.
    -No phisical separation between microlenses is considered.
    - Only those microlenses which are completely inside the main pupil are
    selected for phase recovering
    - When propagation between lenslet and CCD is considered, Fraunhofer
    approximation is considered. Fresnel propagation has as many conditions
    that I have not find a way to make it work with a normal S-H
    configuration. For better accuracy, maybe a ray tracing program will
    show more realistic behaviour. Nevertheless, this software provides
    qualitatively valid results and show a realistic approximation to the
    problem.
    - Only square lenslet and CCD are considered
    - In the case of no propagation between lenslet and CCD, a pad with
    zeros should be implemented to the incoming phase for each microlenses
    in order to avoid edge effects.


Example of use:
	modes=36;
	ZerValues(1:15)=rand(15,1);
	ZerValues(16:modes)=0.1*rand(length(16:modes),1);
	[RecoveredZern]=SHsimulator(36,ZerValues,1000,10e-6,0.780,20,10e-3,1,1e-8,16,1,1,1,1,1,0.1);
%}

%%
close all
if nargin ~= 18 && nargin>0
    error('This function requires 15 inputs, as defined in the documentation');
elseif nargin ==0
    modes=36;ZerValues(1:15)=rand(15,1); ZerValues(16:modes)=0.1*rand(length(16:modes),1);
    resolution=1000; PixelSize=10e-6; LAMBDA=0.78; Lenses=20; focal=10e-3; Prop=1; factor=1e-8; 
    bits=16; MLGeometry=1 ;paint=1;MLCentroidMethod=1; MLSharedArea=1; MLAberration=1; 
    MLfieldDistortion=1;MLvignetting=1; ObjectMagnitude=0; PupilDiameter=4.2; Exposuretime=1e-2; 
    BandWidth=88; QE=0.90; WellCapacity=18e3; PhotonNoise=1; ReadNoise=1; DarkCurrent=0.01; ...
        ReadOutNoise=8;
    fprintf('Function called without values. Using values by defect. Read documentation to include your own values.\n')
end
%Create the configuration of the Shack-Hartmann and microlenses array in a struct.
[SH,ML]=createSH(modes,resolution,PixelSize,LAMBDA,Lenses,focal,Prop,factor,bits,MLGeometry,...
    MLCentroidMethod,MLSharedArea,MLAberration,MLfieldDistortion,MLvignetting,WellCapacity,...
    DarkCurrent,ReadOutNoise,paint);

%Calculate number of photons in reaching detector
nFot = 1000/(10^(ObjectMagnitude/2.5)); %number of photons from the object (fot/s/cm2/Angstrom)
CollectorAreaCm2=pi*((100*PupilDiameter/2)^2); %CollectorArea in cm^2
nFot = round(nFot*Exposuretime*CollectorAreaCm2*BandWidth); 
nFot=nFot*QE;%total photons reaching the detector.

%%
%{
*************************************************************************
********************Calculation of Zernike matrix************************
*************************************************************************
First, search if a Zernike Matrix with the actual number of modes already exist.
%}
if exist ('ZernikeMatrix.mat','file') ~=0 && exist ('config.mat','file') ~=0
    config_=load ('config.mat');
    config=config_.config;
    clear config_
    if config.resolution==SH.resolution && config.modes==SH.modes
        fprintf ('Zernike matrix already exists.\n')
        ZerModo=load ('ZernikeMatrix.mat');
        ZerModo=ZerModo.ZerModo;
    else
        delete('ZernikeMatrix.mat');
        delete('config.mat');
        fprintf ('Configuration has ben modified. Calculating a new Zernike matrix...\n')
        ZerModo=cell(SH.modes,1);
        for i=1:SH.modes
            ZerModo{i}=zernike(i,SH.resolution);
        end
        config=struct('resolution',SH.resolution,'modes',SH.modes); %#ok<NASGU>
        save('ZernikeMatrix.mat', 'ZerModo')
        save('config.mat', 'config')
    end
else
    fprintf('Zernike matrix is not available. Calculating...\n')
    ZerModo=cell(SH.modes,1);
    for i=1:SH.modes
        ZerModo{i}=zernike(i,SH.resolution);
    end
    config=struct('resolution',SH.resolution,'modes',SH.modes);  %#ok<NASGU>
    save('ZernikeMatrix.mat', 'ZerModo')
    save('config.mat', 'config')
end


%%
%{
*************************************************************************
**** Introducing the incoming wavefront in microns (NOLL NOTATION**********
Zernike polyomials and atmospheric turbulence J Op Soc Am. Vol 66, No 3 ,
************************************March 1976****************************
It is suppose to set the 3 first value to zero (piston and tip/tilt)
%}
ZSum=zeros(length(ZerModo{1}));
for i=2:SH.modes
    ZMode=ZerValues(i)*ZerModo{i};
    ZSum=ZSum+ZMode;
end
WF=ZSum.*1e-6;%conversion to meters

%Painting
[SH.pupil4paint]=PaintWavefront(SH,WF);


%%
%{
*************************************************************************
***************Definition of the microlenses array***********************
*************************************************************************
%}
[ML]=MicroLenses(SH,ML,WF);

%{
By defect, Microlenses are perfect shperical lenses according to the
ideal Fraunhofer propagation. If some aberrations has to be added to
microlenses, choose the zernike to be implemented and the amount of microns
in the following form:
ML.Aberration=    %This select the Zernike polynomial to be added in Noll
notation (i.e. ML.Aberration=11)
ML.AberrationAmount=  %Amount of zernike polynomial in meters
(i.e. ML.AberrationAmount= 1e-6)
%}

WFSubpupil=(WF+ML.AmplitudeMask);%Wavefront in each microlent (phase of 
%microlenses itself is not yet taking into account

Lenses=cell(1,length(ML.coor));%Preallocation
for i=1:length(ML.coor)
    Lenses{i}=WFSubpupil(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4));
end



%%
%{
*************************************************************************
*************** Point Spread Function Calculation ***********************
*************************************************************************
%}
PSF=cell(1,length(Lenses));
PSFmaxLocal=zeros(1,length(Lenses));
if ML.Prop==0 %CASE OF NO PROPAGATION BETWEEN MICROLENSES AND CCD. Edge effects will exist, it should be implemented a pad with zeros to avoid it
    fprintf('Propagation between microlenses and CCD is not considered. Microlenses are considered as pure amplitude objects\n')
    for i=1:length(Lenses)
        PF=ML.Pupil.*exp(-1i*SH.k.*Lenses{i});%Pupil function
        PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;%PSF
        PSFmaxLocal(i)=max(max(PSF{i}));
    end
else
    L=length(Lenses{1})*SH.PixelSize; %Size of the region of each microlenses
    for j=1:length(Lenses)
        S = ML.Pupil.*exp(-1i*SH.k.*Lenses{j}); %complex phase screen
        
        if isfield(ML,'Aberration') == 1 %If aberration o each microlens itself is considered
            aberration = zernike(ML.AberrationZernike,size(S,1));
            S = S.*exp(-1i*SH.k*aberration*ML.AberationValue*SH.LAMBDA);
        end
        
        PSF{j} = lensletSimulation(S,L,SH.LAMBDA,ML.focal,SH.PixelSize);
        
        if isfield(ML,'fieldDistortion') == 1
            PSF{j} = applySeidelDistortion(ML.focal,SH.PixelSize,0.1*SH.LAMBDA,size(S,1),PSF{j});
        end
        
        if isfield(ML,'vignetting') == 1
            PSF{j} = applyVignetting(ML.focal,SH.PixelSize,size(S,1),1e-3,PSF{j});
        end
        
        PSFmaxLocal(j)=max(max(PSF{j}));
    end
end



%%
%{
***************************************************************************
****** Quantization of the signal, introduction of photonic noise**********
****************** & introduction of read noise****************************
***************************************************************************
Photon noise, also known as Poisson noise, is a basic form of uncertainty 
associated with the measurement of light, inherent to the quantized nature 
of light and the independence of photon detections. Its expected magnitude 
constitutes the dominant source of image noise except in low-light conditions.
Individual photon detections can be treated as independent events that 
follow a random temporal distribution. As a result, photon counting is a 
classic Poisson process.Photon noise is signal  dependent, and its standard 
deviation grows with the square root of the signal. Contrary to popular 
belief, shot noise experienced by the detector IS related to the QE of the 
detector! Back-illuminated sensors with higher QE yields a better 
Signal/Shot Noise ratio. There is a simple intuitive explanation for this –
 shot noise must be calculated from the signal represented by the number 
of photoelectrons in the sensor (electrons generated from photons falling 
on the sensor), NOT JUST from the number of incoming photons. Therefore, 
if an average of 100 photons hit a pixel, but the sensor has a QE of 50% at
the wavelength of these photons, then an average of 50 photoelectrons will 
be created – the shot noise and Signal/Shot Noise must be calculated from 
this value.
%}
if SH.bits==0
    fprintf('Calculations performed with "double" precision of MATLAB. \n');
    for j=1:length(PSF)
        PSF{j}=PSF{j}/max(max(PSF{j}));
    end
else
    PSFmax=max(max(PSFmaxLocal));
    fprintf('Calculations performed for %2.0f bits. \n',SH.bits);
    for i=1:length(PSF)
        PSF{i}=round((PSF{i}/PSFmax)*((2^SH.bits)-1));
    end
end

% paintingBigPSF.mat included the code that create the PSF "seen" by each 
% microlens, even when there exist energy from surrounding microlenses. It
% also includes the introduction of photonic and read noise.
if ML.radiusPixels*2==length(PSF{1})
    [PSF,idealPSF]=paintingPSF(PSF,ML,SH,nFot,0,0,PhotonNoise,ReadNoise,Exposuretime);
    drawnow();
elseif ML.radiusPixels*2~=length(PSF{1})
    [PSF,idealPSF]=paintingBigPSF(PSF,ML,SH,nFot,0,0,PhotonNoise,ReadNoise,Exposuretime);
    drawnow();
end



%%
%{
**************************************************************************
******************** Centroid calculation*********************************
**************************************************************************
Some problems regarding centroid estimation in Shak-Hartmann:
http://www.ctio.noao.edu/soar/sites/default/files/SAM/archive/5490-123.pdf
Some methods for improve this calculation: "Shack-Hartmann wavefront sensor 
image analysis: a comparison of centroiding methods and image-processing 
techniques"
%}
Xcentroid=zeros(1,length(PSF));
Ycentroid=Xcentroid;
%If you want to visualize each centroid calculation, set the variable
%'display' =1. It is not in the main config because it is usefull only for
%testing
display=0;
if display==1
    figure
end
for i=1:length(PSF)
    [Xcentroid(i),Ycentroid(i)]=CentroidCalculation(PSF{i},ML.CentMethod,display);
end


%Coordinates of the real centroids and the reference centroids are calculated.
ML.RefCentroid=ML.coor+ML.radiusPixels;% Coordinates of the reference centroid.
CentroidLength=zeros(length(ML.coor),2);
CoorCentroid=zeros(length(ML.coor),2);
crosstalks=0;
for i=1:length(ML.RefCentroid)
    CentroidLength(i,1)=Xcentroid(i)-(ML.radiusPixels);
    CentroidLength(i,2)=Ycentroid(i)-(ML.radiusPixels);
    CoorCentroid(i,1)=ML.RefCentroid(i,1)+ CentroidLength(i,1);
    CoorCentroid(i,2)=ML.RefCentroid(i,3)+CentroidLength(i,2);
end

%The following calculates if there exist double spots or crosstalk between
%microlenses.
if isempty(idealPSF)==1
    [IdealXcentroid(i),IdealYcentroid(i)]=CentroidCalculation(idealPSF{i},ML.CentMethod,0);
end
if isempty(idealPSF)==1
    for i=1:length(ML.RefCentroid)
        CentroidLength(i,1)=IdealXcentroid(i)-(ML.radiusPixels);
        CentroidLength(i,2)=IdealYcentroid(i)-(ML.radiusPixels);
        CoorCentroid(i,1)=ML.RefCentroid(i,1)+ CentroidLength(i,1);
        CoorCentroid(i,2)=ML.RefCentroid(i,3)+CentroidLength(i,2);
        if CentroidLength(i,1)>ML.radiusPixels
            crosstalks=crosstalks+1;
        else
            if CentroidLength(i,2)>ML.radiusPixels
                crosstalks=crosstalks+1;
            end
        end
    end
end
if SH.paint==1
    figure
    suptitle('Reference centroid vs calculated centroid')
    subplot(1,2,1)
    set(gcf,'color','w');
    plot(ML.RefCentroid(:,1),ML.RefCentroid(:,3),'o','MarkerSize',2,'MarkerEdgeColor','r')
    hold on
    plot(CoorCentroid(:,1)',CoorCentroid(:,2)','x','MarkerSize',5,'MarkerEdgeColor','b')
    hold on
    for i=1:length(PSF)
        rectangle('Position', [ML.coor(i,1) ML.coor(i,3) ML.radiusPixels*2 ML.radiusPixels*2],'LineWidth', 0.1, 'EdgeColor', 'b');
        %pause(0.01)
    end
    legend('Reference centroid','calculated centroid')
    title('Centroids position')
    xlabel({'pixels';['Number of centroids exceeding its microlent area (crosstalk): ' num2str(crosstalks)]});
    ylabel('pixels')
    set(gca,'ydir','reverse')
    xlim([0 SH.resolution])
    ylim([0 SH.resolution])
    
    subplot(1,2,2)
    quiver(ML.RefCentroid(:,1),ML.RefCentroid(:,3),CentroidLength(:,1),CentroidLength(:,2),0);
    title('Displacement of the real centroid')
    xlabel('pixels')
    ylabel('pixels')
    set(gca,'ydir','reverse')
    xlim([0 SH.resolution])
    ylim([0 SH.resolution])
    drawnow();
end


%%
%{
**************************************************************************
********************** Slope Calculations*********************************
**************************************************************************
%}
delta=double(CentroidLength*SH.PixelSize);
alfax= double(atan(delta(:,1)/ML.focal));
alfay=double(atan(delta(:,2)/ML.focal));
solvingX=ML.erased;
solvingY=ML.erased;

count=1;
for i=1:length(ML.erased)
    if ML.erased(i)==1
        solvingX(i)=alfax(count);
        solvingY(i)=alfay(count);
        count=count+1;
    end
end

if SH.paint==1
    % Building a slopes matrix
    solvingXvec=vec2mat(solvingX',sqrt(length(ML.erased)));
    solvingYvec=vec2mat(solvingY',sqrt(length(ML.erased)));
    solvingX2=solvingXvec;
    solvingY2=solvingYvec;
    for i=length(solvingXvec):-1:1
        for j=length(solvingXvec):-1:1
            if solvingX2(i,j)==0
                solvingX2(i,j)=NaN;
            end
            if solvingY2(i,j)==0
                solvingY2(i,j)=NaN;
            end
        end
    end
    
    figure
    suptitle('Slope or the recovered phase')
    subplot(1,2,1)
    set(gcf,'color','w');
    imshow(solvingX2,[])
    title('x-axis')
    colorbar
    
    subplot(1,2,2)
    imshow(solvingY2,[])
    title('y-axis')
    colorbar
    drawnow();
end

%Joining together in a vector
Solving=zeros(length(solvingX)*2,1);
Solving(1:length(solvingX),1)=solvingX;
Solving(length(solvingX)+1:length(solvingY)*2,1)=solvingY;


%%
%{
**************************************************************************
********************** Slope Calculations*********************************
**************************************************************************
For phase recovering I will use lineal estimation without ligatures,
%using the modal estimation case:"Wavefront Optics for Vision Correction", Guang-ming Dai.


% ****************CALCULATES A VECTOR OF ZERNIKE COEFFICIENTS FOR COMPARATION PURPOSES*****************
% First, it is checked if such vector already exist with the actual configuration
%}
if exist ('SH.mat','file')==2 && exist ('ML.mat','file')==2
    temp=load ('SH.mat');
    temp2=load ('ML.mat');
else
    temp.SH=0;
    temp2.ML=0;
end

if exist ('CalibratedZernike.mat','file')==2 && isequaln(temp.SH,SH)==1 && isequaln(temp2.ML,ML)==1
    fprintf ('Zernike calibration already exists for this configuration.\n')
    CalibratedZernike_=load('CalibratedZernike.mat');
    CalibratedZernike=CalibratedZernike_.CalibratedZernike;
else
    delete('CalibratedZernike.mat');
    delete('SH.mat');
    delete('ML.mat');
    if exist ('CalibratedZernike.mat','file')==2
        fprintf ('Configuration was modified. Calculating new Zernike calibration...\n')
    else
        fprintf('Zernike calibration is not available. Calculating...\n')
    end
    [CalibratedZernike] = ModalZernikeCalibration(SH,ML,ZerModo,length(Solving),nFot);
    save('CalibratedZernike.mat', 'CalibratedZernike');
    save('SH.mat', 'SH');
    save('ML.mat', 'ML');
end


%%
%{
**************************************************************************
************************ Recovering***************************************
**************************************************************************
Only modal recovering has been implanted.
Recovering is made by mean of least square method. It is equivalent to construct the recovering matrix.
%}

RecoveredZern = lsqr((CalibratedZernike(:,1:SH.modes)*1e-6)/SH.factor,Solving,1e-10,500);
%Recovered wavefront
RecWF=zeros(length(ZerModo{1}));
for i=4:SH.modes
    ZModeRec=RecoveredZern(i)*ZerModo{i};
    RecWF=RecWF+ZModeRec;
end
RecWF=RecWF.*1e-6;%conversion to meters.
paintWFs(WF,RecWF,SH.pupil4paint);

%%
%{
*************************************************************************
***********RESIDUAL ANALYSIS********************************************
**************************************************************************
 *************************************************************************
A good book for residual theory is the thesis dissertation of Justo Arines
%}
Residual=SH.pupil.*(WF-RecWF);

%Pupil function
PFOrig = SH.pupil.*exp(sqrt(-1)*SH.k.*WF);
PFRec = SH.pupil.*exp(sqrt(-1)*SH.k*RecWF);
PFDif=SH.pupil.*exp(sqrt(-1)*SH.k*(RecWF*0));
PFRes=SH.pupil.*exp(sqrt(-1)*SH.k*Residual);


%PSF
PSFOrig=fft2(PFOrig); %Two-dimensional discrete Fourier Transform.
PSFOrig=fftshift(PSFOrig);
PSFOrig=PSFOrig.*conj(PSFOrig);

PSFRec=fft2(PFRec); %Two-dimensional discrete Fourier Transform.
PSFRec=fftshift(PSFRec);
PSFRec=PSFRec.*conj(PSFRec);

PSFDif=fft2(PFDif); %Two-dimensional discrete Fourier Transform.
PSFDif=fftshift(PSFDif);
PSFDif=PSFDif.*conj(PSFDif);

PSFRes=fft2(PFRes); %Two-dimensional discrete Fourier Transform.
PSFRes=fftshift(PSFRes);
PSFRes=PSFRes.*conj(PSFRes);


%OTF & MTF
OTFOrig=fft2(PSFOrig);%OTF
MTFOrig = abs(OTFOrig);
MTFOrig=fftshift(MTFOrig);
MTFOrig = MTFOrig./max(max(MTFOrig));

OTFRec=fft2(PSFRec);
MTFRec = abs(OTFRec);
MTFRec=fftshift(MTFRec);
MTFRec = MTFRec./max(max(MTFRec));

OTFRes=fft2(PSFRes);
MTFRes = abs(OTFRes);
MTFRes=fftshift(MTFRes);
MTFRes = MTFRes./max(max(MTFRes));


if SH.paint==1
    paintPSFs(PSFOrig,PSFRec,PSFRes)
    PaintConvolution(PSFOrig,PSFRec,PSFRes)
    paintMTF(SH,MTFOrig,MTFRec,MTFRes)
end


% *************************************************************************
% ******************* Quality metrics**************************************
% *************************************************************************

StRes=max(max(PSFRes))/max(max(PSFDif));
fprintf('Strhel ratio of residual = %2.5f (ideal =>0.8)\n',StRes);

StOrig=max(max(PSFOrig))/max(max(PSFDif));
StRec=max(max(PSFRec))/max(max(PSFDif));
fprintf('Difference between Strhel ratio of the incoming phase and the recovered one is %2.5f\n',StOrig-StRec);

RMSOrig=sqrt(sum(ZerValues(4:length(ZerValues)).^2));
RMSRec=sqrt(sum(RecoveredZern(4:length(RecoveredZern)).^2));
RecoveredZern=RecoveredZern';
RMSDifference=sqrt(sum(abs(ZerValues(4:length(ZerValues))-RecoveredZern(4:length(RecoveredZern))).^2));
fprintf('Difference between the RMS of the incoming phase (%2.5f) and the recovered one (%2.5f) is = %2.5f (%2.2f of error in percentage).\n',RMSOrig,RMSRec,RMSDifference,100-(RMSRec*100/RMSOrig));


figure
plot(4:SH.modes,ZerValues(4:SH.modes))
hold on
plot(4:SH.modes,RecoveredZern(4:SH.modes))
set(gcf,'color','w');
legend('Original','Recovered');
xlabel('Zernike mode')
ylabel('Value in microns')
title('Original Zernike coefficient Vs Recovered')
drawnow();

