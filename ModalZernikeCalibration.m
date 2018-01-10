function [CalibratedZernike] = ModalZernikeCalibration(SH,ML,ZerModo,large,nFot)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function characterize the Shack-Hartmann by mean of calculating the output of each individual Zernike mode through the system.

OUTPUTS:
CalibratedZernike: A matrix where each column is a vector with the result of use each individual zernike mode through the system.
%}
CalibratedZernike=zeros(large,SH.modes);
for w=1:SH.modes;
    ZMODE=ZerModo{w};
    W=ZMODE.*SH.factor;%conversion to meters
    WFSubpupil=W.*ML.AmplitudeMask;% Wavefront seen through microlenses array
    Lens=cell(1,length(ML.coor));%Preallocation
    PSF=cell(1,length(Lens));%Preallocation
    PSFmaxLocal=zeros(1,length(Lens));
    
    for i=1:length(ML.coor)
        Lens{i}=WFSubpupil(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4));
    end
    
    
    %CALCULATING THE PSF OF EACH MICROLENSE
    if ML.Prop==0 % CASE OF NO PROPAGATION BETWEEN MICROLENSES AND CCD
        for i=1:length(Lens)
            PF=ML.Pupil.*exp(-1i*SH.k.*Lens{i});
            PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;
            PSFmaxLocal(i)=max(max(PSF{i}));
        end
    else %PROPAGATION BETWEEN MICROLENSES AND CCD EXISTS
        L=length(Lens{1})*SH.PixelSize; %Size of the region of each microlenses
        for j=1:length(Lens)
            S = ML.Pupil.*exp(-1i*SH.k.*Lens{j}); %complex phase screen
        
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
    
    for j=1:length(PSF)
        PSF{j}=PSF{j}/max(max(PSF{j}));
    end
    
    %For calibration no quantization or noise is taken into account.
    if ML.radiusPixels*2==length(PSF{1})
        [PSF,~]=paintingPSF(PSF,ML,SH,nFot,1,1,0,0,0);
    elseif ML.radiusPixels*2~=length(PSF{1})
        [PSF,~]=paintingBigPSF(PSF,ML,SH,nFot,1,1,0,0,0);
    end
    
    
    
    Xcentroid=zeros(1,length(PSF));
    Ycentroid=Xcentroid;
    display=0;
    for i=1:length(PSF)
        [Xcentroid(i),Ycentroid(i)]=CentroidCalculation(PSF{i},ML.CentMethod,display);
    end
    
    %Coordinates of the real centroids and the reference centroids are calculated.
    ML.RefCentroid=ML.coor+ML.radiusPixels;% Coordinates of the reference centroid.
    ML.CentroidLength=zeros(length(ML.coor),2);
    CoorCentroid=zeros(length(ML.coor),2);
    for i=1:length(ML.RefCentroid)
        ML.CentroidLength(i,1)=Xcentroid(i)-(ML.radiusPixels);
        ML.CentroidLength(i,2)=Ycentroid(i)-(ML.radiusPixels);
        CoorCentroid(i,1)=ML.RefCentroid(i,1)+ ML.CentroidLength(i,1);
        CoorCentroid(i,2)=ML.RefCentroid(i,3)+ML.CentroidLength(i,2);
    end
    
        
    delta=double(ML.CentroidLength*SH.PixelSize);
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
    
    
    %Joining together in a vector
    Solving=zeros(length(solvingX)*2,1);
    Solving(1:length(solvingX),1)=solvingX;
    Solving(length(solvingX)+1:length(solvingY)*2,1)=solvingY;
    CalibratedZernike(:,w)=Solving;
end



