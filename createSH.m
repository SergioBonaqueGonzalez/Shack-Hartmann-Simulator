function [SH,ML]=createSH(modes,resolution,PixelSize,LAMBDA,Lenses,...
    focalML,Prop,factor,bits,MLGeometry,MLCentroidMethod,MLSharedArea,...
    MLAberration,MLfieldDistortion,MLvignetting,WellCapacity,DarkCurrent,ReadOutNoise,paint)

SH.modes=modes; 
SH.resolution=resolution; 
SH.PixelSize=PixelSize; 
SH.LAMBDA=LAMBDA*1e-6; 
SH.factor=factor; 
SH.bits=bits; 
SH.paint=paint; 
SH.k=2*pi/SH.LAMBDA; 
SH.SensorSize=SH.PixelSize*SH.resolution;
SH.radius=SH.SensorSize/2;
[SH.pupil,SH.rho]=PupilCreator(SH.resolution); %We are using circular pupils.
SH.WellCapacity=WellCapacity;
SH.DarkCurrent=DarkCurrent;
SH.ReadOutNoise=ReadOutNoise;

ML.focal=focalML; 
ML.Prop=Prop; 
ML.Lenses=Lenses; 
ML.spacing=floor(SH.resolution/ML.Lenses(1)); %Asigned pixels in CCD to each microlense
surplus_=SH.resolution-ML.Lenses*ML.spacing;
ML.surplus=[floor(surplus_) ceil(surplus_)];
ML.Geometry=MLGeometry;
ML.radiusPixels=ML.spacing/2; %Radius of microlenses in pixels
ML.FresnelNumber=((SH.SensorSize/2)^2)/(SH.LAMBDA*ML.focal); %Fresnel expression describes diffraction under the paraxial assumption, where only rays that make a small angle (< ~0.1 rad) relative to the optical axis are considered.
ML.CentMethod=MLCentroidMethod; %Used method for centroid calculations
ML.SharedArea=MLSharedArea;
ML.fieldDistortion=MLfieldDistortion;
ML.vignetting=MLvignetting;

ML.Aberration=MLAberration; %Flag that indicates if microlenses has its own aberrations.
% Default values for aberrations of the microlenses
if ML.Aberration==1
    ML.AberrationZernike=11; %Zernike coefficient
    ML.AberationValue=0.1; %Value in microns
end

if  rem(surplus_,2)~=0
    fprintf('WARNING:this configuration does not allow a symmetrical disposition of microlenses. \nYou can try to ajust the resolution to an odd number and to select another number of microlenses.\n')
end
