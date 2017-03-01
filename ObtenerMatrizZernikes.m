function [VectorZernikes] = ObtenerMatrizZernikes(TamPixel,focalML,Propagacion,lambda,k,MicroLentesMascara,radioMLpxs,CoorBuenas,PupilaML,NFresnel,CentroideReferencia,Escala,modos,ZerModo,Eliminadas,longitud,factor,bits)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
%This function characterize the Shack-Hartmann by mean of calculating the output of each individual Zernike mode through the system.

%INPUTS:
%All inputs are self-called from the main function.
%TamPixel: size of pixel in the CCD.
%focalML: focal length of microlenses.
%Propagacion: flag which indicates if propagation between microlenses and CCD should be considered or not. 
%lambda: wavelength in meters.
%k= 2*pi/lambda
%MicroLentesMascara= mask of the microlenses
%radioMLpxs= radius of each microlense in pixels.
%CoorBuenas= Coordinates of the upper left corner of each valid microlense. A valid microlense is the one which is completely inside the whole pupil of the system.
%PupilaML= tiny binary pupil of each microlense
%NFresnel: Fresnel number of microlenses
%CentroideReferencia: coordinates of the refence centroids
%Escala: Each microlense is virtually resized to its double for PSF calculation. It can be changed by mean of this variable.
%modos: number of Zernike modes to be characterized
%ZerModo: Matrix which contain each Zernike mode in the proper size
%Eliminadas:Coordinates of the upper left corner of each invalid microlense. For example, discarded microlenses because they are incomplete.
%longitud: size of the expected vector size for each mode. For pre-allocation purposes.
%factor:is the value with which Zernike modes are characterized
%bits: number of bits of the CCD

%OUTPUTS:
%VectorZernikes: A matrix where each column is a vector with the result of use each individual zernike mode through the system.

VectorZernikes=zeros(longitud,modos);
PSFmax=zeros(1,modos);
parfor w=1:modos;
    ZMODO=ZerModo{w};
    W=ZMODO.*factor;%conversion to meters
    WFSubpupil=W.*MicroLentesMascara;% Wavefront seen through microlenses array 
    Lente=cell(1,length(CoorBuenas));%Preallocation
    PadLente=cell(1,length(CoorBuenas));%Preallocation
    
    for i=1:length(CoorBuenas)
        Lente{i}=WFSubpupil(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4));% Spliting the wavefront in the microlenses size
        PadLente{i}=zeros(radioMLpxs*Escala);% Resize the matrix 
        PadLente{i}(radioMLpxs*2+1:end-radioMLpxs*2,radioMLpxs*2+1:end-radioMLpxs*2)= Lente{i};
    end

    %CALCULATING THE PSF OF EACH MICROLENSE
    if Propagacion==0 % CASE OF NO PROPAGATION BETWEEN MICROLENSES AND CCD
        PSF=cell(1,length(PadLente));%Prealoco los arrays
        PSFmaxLocal=zeros(1,length(PadLente));
        for i=1:length(PadLente)
            PF=PupilaML.*exp(sqrt(-1)*k.*PadLente{i});
            PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;
            PSFmaxLocal(i)=max(max(PSF{i}));
        end
    else %PROPAGATION BETWEEN MICROLENSES AND CCD EXISTS
        [~,a]=size(PadLente{1});
        L=a*TamPixel;%Size of the region of each microlenses
        if NFresnel> 0.5 && NFresnel<1;
            PSF=cell(1,length(PadLente));
            PSFmaxLocal=zeros(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensity of image
                PSFmaxLocal(i)=max(max(PSF{i}));
            end
        elseif NFresnel>= 1
            PSF=cell(1,length(PadLente));
            PSFmaxLocal=zeros(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensity of image
                PSFmaxLocal(i)=max(max(PSF{i}));
            end
        else
            PSF=cell(1,length(PadLente));
            PSFmaxLocal=zeros(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFraunhofer(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensity of image
                PSFmaxLocal(i)=max(max(PSF{i}));
            end
        end
    end
    
    PSFmax(w)=max(max(PSFmaxLocal));
end

PSFmaxTotal=max(max(PSFmax));
%First I calculate the total maximum value. After that I repeated the process, scaling to this value. 

parfor w=1:modos;
    ZMODO=ZerModo{w};
    W=ZMODO.*factor;
    WFSubpupil=W.*MicroLentesMascara;
    Lente=cell(1,length(CoorBuenas));
    PadLente=cell(1,length(CoorBuenas));
    for i=1:length(CoorBuenas)
        Lente{i}=WFSubpupil(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4))
        PadLente{i}=zeros(radioMLpxs*Escala);
        PadLente{i}(radioMLpxs*2+1:end-radioMLpxs*2,radioMLpxs*2+1:end-radioMLpxs*2)= Lente{i};
    end

    %CALCULATING THE PSF OF EACH MICROLENSE
    if Propagacion==0 % CASE OF NO PROPAGATION BETWEEN MICROLENSES AND CCD
        PSF=cell(1,length(PadLente));%Prealoco los arrays
        for i=1:length(PadLente)
            PF=PupilaML.*exp(sqrt(-1)*k.*PadLente{i});%Pupil function
            PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;%PSF
        end
    else %PROPAGATION BETWEEN MICROLENSES AND CCD EXISTS
        [~,a]=size(PadLente{1});
        L=a*TamPixel;%Size of the region of each microlenses
        if NFresnel> 0.5 && NFresnel<1;
            PSF=cell(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensity of image
            end
        elseif NFresnel>= 1
            PSF=cell(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensity of image
            end
        else
            PSF=cell(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFraunhofer(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensity of image
            end
        end
    end
    
    centroideRealX=zeros(1,length(PSF));
    centroideRealY=zeros(1,length(PSF));
    
    for i=1:length(PSF)
        [centroideRealX(i),centroideRealY(i)]=CalCentroides(PSF{i},1,0);
    end
    
    LongitudCentroide=zeros(length(CoorBuenas),2);
    CoorCentroide=zeros(length(CoorBuenas),2);
    for i=1:length(CentroideReferencia)
        LongitudCentroide(i,1)=centroideRealX(i)-(Escala*radioMLpxs/2);
        LongitudCentroide(i,2)=centroideRealY(i)-(Escala*radioMLpxs/2);
        CoorCentroide(i,1)=CentroideReferencia(i,1)+ LongitudCentroide(i,1);%coordinates of the centroid
        CoorCentroide(i,2)=CentroideReferencia(i,3)+LongitudCentroide(i,2);
    end
    
    deltas=double(LongitudCentroide*TamPixel); 
    alfax= double(atan(deltas(:,1)/focalML));
    alfay=double(atan(deltas(:,2)/focalML));
    
    %each alfa is placed in its correct place in the initial matrix. 
    solucionx=Eliminadas;
    soluciony=Eliminadas;
    contador=1;
    for i=1:length(Eliminadas)
        if Eliminadas(i)==1
            solucionx(i)=alfax(contador);
            soluciony(i)=alfay(contador);
            contador=contador+1;
        end
    end
    
    %Results are jointed in a vector: 
    solucionbi=zeros(length(solucionx)*2,1);
    solucionbi(1:length(solucionx),1)=solucionx;
    solucionbi(length(solucionx)+1:length(soluciony)*2,1)=soluciony;
    VectorZernikes(:,w)=solucionbi;
end

