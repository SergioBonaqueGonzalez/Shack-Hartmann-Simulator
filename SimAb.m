function [ZerRecuperados]=SimAb(modos,c,resolucion,TamPixel,LAMBDA,NLentes,focalML,Propagacion,factor,bits,pintar)

%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This program simulates a Shack-Hartmann sensor. 

%INPUTS
%***modos= number of Zernike modes inteded to be recovered (ej. modos=10)
%***c = vector that contains the Zernike coefficients of the incoming phase in Noll notation. (ej. c=rand(10,1)
%***resolucion=resolution of CCD (Only square and odd CCDs are suported) (ej. 1024)
%***TamPixel= Pixel size of the CCD (ej. 1.471e-6)
%***LAMBDA= wavelength of the simulations in microns (ej. 0.780)
%***NLentes= number of microlenses in a row in the microlenses array (ej.41)
%***focalML= Focal length of the microlenses in meters.
%***Propagacion= flag that indicates if propagation should be considered between microlenses and CCD =0 NO propagation; =1 Propagation 
%***factor= the value with which zernize coefficients are characterized through the system. ej=1e-8
%***bits= number of bits of the CCD 
%***pintar= dummy flag which indicates if all figures should be painted 1=YES, 2=NO

%%Example of use:
%modos=36;
%c(1:15)=rand(15,1);
%c(16:modos)=0.1*rand(length(16:modos),1);
%[ZerRecuperados]=SimAb(modos,c,2025,1.471e-6,0.780,41,1e-3,0,1e-8,16,0);


tic



% *************************************************************************
% ************************SENSOR FEATURES *********************************
% *************************************************************************
% *************************************************************************
% close all
% clear all
% clc
c=c(:);%This line avoid annoying messages when painting
warning('off', 'Images:initSize:adjustingMag'); 


if  rem(resolucion,2)==0
    error('Error. Resolution of the CCD should be odd.')
end

espaciado=floor(resolucion/NLentes(1)); %Asigned pixels in CCD to each microlense
Resto=resolucion-NLentes(1)*espaciado;

if  rem(Resto,2)~=0
    error('Error.This configuration does not allow a symmetric configuration of microlenses array. Please, change the number of microlenses or the resolution.')
end



lambda = LAMBDA*1e-6;%Wavelength in meters
k = 2*pi/lambda;
TamanoSensor=TamPixel*resolucion;
SeparacionML=0; % Distance between microlenses. Now it is suppossed to be zero. This is something to be implemented in the future

Radio=TamanoSensor/2;%Radius of the pupil
[pupil,rho]=CrearPupila(resolucion);


% *************************************************************************
% ************************Calculation Zernike matrix**********************
% *************************************************************************
% *************************************************************************
%First, search for the Zernike Matrix with the actual configuration of the simulator.

if exist ('MatrizZernikes.mat','file') ~=0 && exist ('config.mat','file') ~=0
    est=load ('config.mat');
    config=est.config;
    clear est
    if config.resolucion==resolucion && config.modos==modos 
        fprintf ('Zernike matrix already exists.\n')
        ZerModo=load ('MatrizZernikes.mat');
        ZerModo=ZerModo.ZerModo;
    else
        delete('MatrizZernikes.mat');
        delete('config.mat');
        fprintf ('Configuration has ben modified. Calculating a new Zernike matrix...\n')
        ZerModo=cell(modos,1);
        for i=1:modos
            ZerModo{i}=zernike(i,resolucion);
        end
        config=struct('resolucion',resolucion,'modos',modos); %#ok<NASGU>
        save('MatrizZernikes.mat', 'ZerModo')
        save('config.mat', 'config')
    end
else
    fprintf('Zernike matrix not available. Calculating...\n')
    ZerModo=cell(modos,1);
    for i=1:modos
        ZerModo{i}=zernike(i,resolucion);
    end
    config=struct('resolucion',resolucion,'modos',modos);  %#ok<NASGU>
    save('MatrizZernikes.mat', 'ZerModo')
    save('config.mat', 'config')
end

clear config



% *************************************************************************
% *********** Introducing the incoming wavefront in microns (NOLL NOTATION**********
%Zernike polyomials and atmospheric turbulence J Op Soc Am. Vol 66, No 3 , March 1976**********************
% *************************************************************************
ZSUMA=zeros(length(ZerModo{1}));

%It is suppose to set the 3 first value to zero.
for i=2:modos
    ZMODO=c(i)*ZerModo{i};
    ZSUMA=ZSUMA+ZMODO;
end
W=ZSUMA.*1e-6;%conversion to meters


                                              

                                            %Painting
                                            [pupilpintar]=PintarFrenteOnda(rho,W,pintar);
                                            
                                            
% *************************************************************************
% ***************DEFINITION OF THE MICROLENSES ARRAY***********************
%**************************************************************************
% *************************************************************************          
[MicroLentesMascara,radioMLpxs,CoorBuenas,Eliminadas]=MicroL(NLentes,resolucion,SeparacionML,pupil,pupilpintar,W,pintar);


% *************************************************************************
% ***********SAMPLING GOODNESS IN EACH MICROLENT**************************
%**************************************************************************
% *************************************************************************     
RadioMLm=TamPixel*radioMLpxs;%Radius in mmm
FML=focalML/(RadioMLm*2);%F number of each microlent
Deltax=FML*lambda; 
MuestreoNecesario=RadioMLm*2/Deltax;

if MuestreoNecesario<=radioMLpxs*2
    fprintf ('Sampling of microlenses is good enough.\n')
else
    fprintf ('ERROR: sampling of microlenses is inadequate.\n')
end

Lente=cell(1,length(CoorBuenas));%Preallocation
PadLente=cell(1,length(CoorBuenas));%Preallocation


% *************************************************************************
% ***********CALCULATING IMAGE OF EACH MICROLENT**********************
%**************************************************************************
% *************************************************************************     
WFSubpupil=W.*MicroLentesMascara;%Wavefront through each microlent
%Each microlent is managed separately
Escala=6; % This value imply an extra space around each microlent when calculating the PSF in order to avoid borders effects. The optimum value is the double of the original.
for i=1:length(CoorBuenas)
        Lente{i}=WFSubpupil(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4));
        PadLente{i}=zeros(radioMLpxs*Escala);
        PadLente{i}(radioMLpxs*2+1:end-radioMLpxs*2,radioMLpxs*2+1:end-radioMLpxs*2)= Lente{i};
end

%Defining the pupil of each microlent in binary.
PupilaML=PadLente{1};
[a,b]=size(PupilaML);
for i=1:a
    for j=1:b
        if PupilaML(i,j)~=0
           PupilaML(i,j)=1;
        end
    end
end

%CALCULATING THE PSF OF EACH MICROLENT
NFresnel=((TamanoSensor/2)^2)/(lambda*focalML); %Fresnel expression describes diffraction under the paraxial assumption, where only rays that make a small angle (< ~0.1 rad) relative to the optical axis are considered. 
if Propagacion==0 %CASE OF NO PROPAGATION BETWEEN MICROLENSES AND CCD
    fprintf('No se considera propagacion entre microlentes y CCD.\n')
    PSF=cell(1,length(PadLente));%Prealoco los arrays
    PSFmaxLocal=zeros(1,length(PadLente));
    for i=1:length(PadLente)
        PF=PupilaML.*exp(sqrt(-1)*k.*PadLente{i});%Pupil function
        PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;%PSF
        PSFmaxLocal(i)=max(max(PSF{i}));
    end
else 
    [~,a]=size(PadLente{1});
    L=a*TamPixel;
    if NFresnel> 0.5 && NFresnel<1;
        fprintf('Numero de Fresnel =%2.0f.\n',NFresnel);
        fprintf('Aware, it is in the border of Fresnel region.\n')
        PSF=cell(1,length(PadLente));
        PSFmaxLocal=zeros(1,length(PadLente));
        for j=1:length(PadLente)
            S = exp(k*1i.*PadLente{j}); %complex phase screen
            propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
            PSF{j} = (abs(propagada).^2);
            PSFmaxLocal(j)=max(max(PSF{j}));
        end
    elseif NFresnel>= 1
        fprintf('Numero de Fresnel =%2.0f.\n',NFresnel);
        fprintf ('Propagation in Fresnel regime.\n')
        PSF=cell(1,length(PadLente));
        PSFmaxLocal=zeros(1,length(PadLente));
        for j=1:length(PadLente)
            S = exp(k*1i.*PadLente{j}); %complex phase screen
            propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
            PSF{j} = (abs(propagada).^2);
            PSFmaxLocal(j)=max(max(PSF{j}));
        end
    else
        fprintf('Numero de Fresnel =%2.0f.\n',NFresnel);
        fprintf ('Propagation in Fraunhofer regime.\n')
        PSF=cell(1,length(PadLente));
        PSFmaxLocal=zeros(1,length(PadLente));
        for j=1:length(PadLente)
            S = exp(k*1i.*PadLente{j}); %complex phase screen
            propagada=propagacionFraunhofer(S.*PupilaML,L,lambda,focalML);
            PSF{j} = (abs(propagada).^2);
            PSFmaxLocal(j)=max(max(PSF{j}));
        end
    end
end

%Quantization according the bits of the CCD
if bits==0
    fprintf('Calculations performed with "double" precision of MATLAB. \n');
    for j=1:length(PSF)
        PSF{j}=PSF{j}/max(max(PSF{j}));
    end
else
    
    PSFmax=max(max(PSFmaxLocal));
    fprintf('Calculations performed for %2.0f bits. \n',bits);
    for i=1:length(PSF)
        PSF{i}=round((PSF{i}/PSFmax)*((2^bits)-1));
    end
end


                                            %Ppainting
                                            if pintar==1
                                                pintarPSFMatriz(PSF,MicroLentesMascara,radioMLpxs,CoorBuenas,pupilpintar,Escala)
                                                drawnow();
                                            end
            
% *************************************************************************
% ***********CALCULO DEL CENTROIDE*****************************************
%**************************************************************************
% *************************************************************************   
%This is a very simple simulation which does not include noise. 
%Some problems regarding centroid estimation in Shak-Hartmann:http://www.ctio.noao.edu/soar/sites/default/files/SAM/archive/5490-123.pdf
%Some methods for improve this calculataion: "Shack-Hartmann wavefront sensor image analysis: a comparison of centroiding methods and image-processing techniques"

centroideRealX=zeros(1,length(PSF));
centroideRealY=zeros(1,length(PSF));

%CalCentroides function calculates the centroid coordinates with sub-pixel precision. 
%In real phase sensor with microlenses arrays, there exists interferences between microlenses apertures, creating one or more secondary spots, finally affecting the calculation. The error term associate to this fact is  =(theta(t)-theta'(t))/theta(t); being tetha the
% the real tilt angle and theta' the calculated oneel calculado ("Analysis on ShackHartmann
%wave-front sensor with Fourier optics"). However, in this case each microlent produces its image separately, so this issue is not taken into account. 
display=0; %1= centroids are painted, 0=no 
if display==1
    figure
end
for i=1:length(PSF)
    [centroideRealX(i),centroideRealY(i)]=CalCentroides(PSF{i},1,display);
end


%Coordinates of the real centroids and the reference centroids are calculated. 
CentroideReferencia=CoorBuenas+radioMLpxs;% Coordinates of the reference centroid.
errores=0;
LongitudCentroide=zeros(length(CoorBuenas),2);
CoorCentroide=zeros(length(CoorBuenas),2);
for i=1:length(CentroideReferencia)
    LongitudCentroide(i,1)=centroideRealX(i)-(Escala*radioMLpxs/2);
    LongitudCentroide(i,2)=centroideRealY(i)-(Escala*radioMLpxs/2);
    CoorCentroide(i,1)=CentroideReferencia(i,1)+ LongitudCentroide(i,1);
    CoorCentroide(i,2)=CentroideReferencia(i,3)+LongitudCentroide(i,2);
    %This introduces a count in order to know where double spots are produced
    if LongitudCentroide(i,1)>radioMLpxs
        errores=errores+1;
    else
        if LongitudCentroide(i,2)>radioMLpxs
            errores=errores+1;
        end
    end
end


                                            %Painting
                                            if pintar==1
                                                figure
                                                suptitle('Reference centroid vs calculated centroid')
                                                subplot(1,2,1)
                                                set(gcf,'color','w');
                                                plot(CentroideReferencia(:,1),CentroideReferencia(:,3),'o','MarkerSize',2,'MarkerEdgeColor','r')
                                                hold on
                                                plot(CoorCentroide(:,1)',CoorCentroide(:,2)','x','MarkerSize',5,'MarkerEdgeColor','b')
                                                hold on
                                                for i=1:length(PSF)
                                                    rectangle('Position', [CoorBuenas(i,1) CoorBuenas(i,3) radioMLpxs*2 radioMLpxs*2],'LineWidth', 0.1, 'EdgeColor', 'b');
                                                    %pause(0.01)
                                                end
                                                legend('Reference centroid','calculated centroid')
                                                title('Centroids position') 
                                                xlabel({'pixels';['Number of centroids exceeding its microlent area (crosstalk): ' num2str(errores)]});
                                                ylabel('pixels')
                                                set(gca,'ydir','reverse')
                                                xlim([0 resolucion])
                                                ylim([0 resolucion])
                                                
                                                subplot(1,2,2)
                                                quiver(CentroideReferencia(:,1),CentroideReferencia(:,3),LongitudCentroide(:,1),LongitudCentroide(:,2),0);
                                                title('Displacement of the real centroid') 
                                                xlabel('pixels')
                                                ylabel('pixels')
                                                set(gca,'ydir','reverse')
                                                xlim([0 resolucion])
                                                ylim([0 resolucion])
                                                drawnow();
                                            end

%If a centroid exceed the area of the CCD "asigned" to its microlent, it is a problem when calculating. Here, we are simulating an ideal Shack-Hartmann, able to find out where each PSF comes from.


% *************************************************************************
% ***********SLOPE CALCULATION***************************************
%**************************************************************************
% *************************************************************************  
deltas=double(LongitudCentroide*TamPixel); 

alfax= double(atan(deltas(:,1)/focalML));
alfay=double(atan(deltas(:,2)/focalML));

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

      
                                             %Painting
                                            if pintar==1
                                                % Building a slopes matrix
                                                solucionxvec=vec2mat(solucionx',sqrt(length(Eliminadas)));
                                                solucionyvec=vec2mat(soluciony',sqrt(length(Eliminadas)));
                                                solucionx2=solucionxvec;
                                                soluciony2=solucionyvec;
                                                for i=length(solucionxvec):-1:1
                                                    for j=length(solucionxvec):-1:1
                                                        if solucionx2(i,j)==0
                                                           solucionx2(i,j)=NaN;
                                                        end
                                                        if soluciony2(i,j)==0
                                                           soluciony2(i,j)=NaN;
                                                        end
                                                    end
                                                end

                                                figure
                                                suptitle('Slope or the recovered phase')
                                                subplot(1,2,1)
                                                set(gcf,'color','w');
                                                imshow(solucionx2,[])
                                                title('x-axis') 
                                                colorbar

                                                subplot(1,2,2)
                                                imshow(soluciony2,[])
                                                title('y-axis') 
                                                colorbar
                                                drawnow();

                                                clear solucionx2
                                                clear soluciony2
                                            end

%Joining together in a vector
solucionbi=zeros(length(solucionx)*2,1);
solucionbi(1:length(solucionx),1)=solucionx;
solucionbi(length(solucionx)+1:length(soluciony)*2,1)=soluciony;

                                            

% *************************************************************************
% ***********ZERNIKE RECOVERING**********************************
%**************************************************************************
% *************************************************************************  
%For phase recovering I will use lineal estimation without ligatures,
%using the modal estimation case:"Wavefront Optics for Vision Correction", Guang-ming Dai.

% ****************CALCULATES A VECTOR OF ZERNIKE COEFFICIENTS FOR COMPARATION PURPOSES*****************
% First, it is checked if such vector already exist with the actual configuration

if exist ('VectorZernikes.mat','file') ~=0
    est=load ('VectorZernikes.mat');
    configpre=est.configuracion;
    clear est
    if configpre.resolucion==resolucion && configpre.TamPixel==TamPixel && configpre.focalML == focalML && configpre.Propagacion==Propagacion && configpre.lambda==lambda && configpre.NLentes == NLentes(1) && configpre.modos==modos && configpre.Escala==Escala && configpre.factor==factor && configpre.bits==bits
        fprintf ('Zernike vector already exists for this configuration.\n')
        VectorZernikes=configpre.VectorZernikes;
    else
        delete('VectorZernikes.mat');
        fprintf ('Configuration was modified. Calculating new vector of Zernikes...\n')
        [VectorZernikes] = ObtenerMatrizZernikes(TamPixel,focalML,Propagacion,lambda,k,MicroLentesMascara,radioMLpxs,CoorBuenas,PupilaML,NFresnel,CentroideReferencia,Escala,modos,ZerModo,Eliminadas,length(solucionbi),factor,bits);
        configuracion = struct('resolucion',resolucion,'TamPixel',TamPixel,'focalML',focalML,'Propagacion',Propagacion,'lambda',lambda,'NLentes',NLentes(1),'modos',modos,'VectorZernikes',VectorZernikes,'Escala',Escala,'factor',factor,'bits',bits);%#ok<NASGU>
        save('VectorZernikes.mat', 'configuracion')
    end
else
    fprintf('Zernike vector is not available. Calculating new Zernike vector...\n')
    [VectorZernikes] = ObtenerMatrizZernikes(TamPixel,focalML,Propagacion,lambda,k,MicroLentesMascara,radioMLpxs,CoorBuenas,PupilaML,NFresnel,CentroideReferencia,Escala,modos,ZerModo,Eliminadas,length(solucionbi),factor,bits);
    configuracion = struct('resolucion',resolucion,'TamPixel',TamPixel,'focalML',focalML,'Propagacion',Propagacion,'lambda',lambda,'NLentes',NLentes(1),'modos',modos,'VectorZernikes',VectorZernikes,'Escala',Escala,'factor',factor,'bits',bits);%#ok<NASGU>                                 
    save('VectorZernikes.mat', 'configuracion')
end

clear configuracion

%Recovering is made by mean of least square method. It is equivalent to construct the recovering matrix.
ZerRecuperados = lsqr((VectorZernikes(:,1:modos)*1e-6)/factor,solucionbi,1e-10,500);


%Recovered wavefront
ZSUMARec=zeros(length(ZerModo{1}));
for i=4:modos
    ZMODORec=ZerRecuperados(i)*ZerModo{i};
    ZSUMARec=ZSUMARec+ZMODORec;
end


                                                pintarWFs(W,ZSUMARec.*1e-6,pupilpintar);
                                                
%Pupil is enlarged when calculating PSF for avoid borders effects
ZTOTALRec=zeros(resolucion*2+1);
ZTOTALPad=ZTOTALRec;
pupilPad=ZTOTALRec;
ZTOTALRec((resolucion/2)+1.5:end-(resolucion/2)-0.5,(resolucion/2)+1.5:end-(resolucion/2)-0.5)=ZSUMARec;
ZTOTALPad((resolucion/2)+1.5:end-(resolucion/2)-0.5,(resolucion/2)+1.5:end-(resolucion/2)-0.5)=ZSUMA;
pupilPad((resolucion/2)+1.5:end-(resolucion/2)-0.5,(resolucion/2)+1.5:end-(resolucion/2)-0.5)=pupil;


WRec=ZTOTALRec.*1e-6;%conversion to meters.
WPad=ZTOTALPad.*1e-6;
                      

% *************************************************************************
% ***********RESIDUAL ANALYSIS********************************************
%**************************************************************************
% *************************************************************************                                            
%A good book for residual theory is the thesis dissertation of Justo Arines
Residuo=pupilPad.*(WPad-WRec);

%Pupil function
PFOrig = pupilPad.*exp(sqrt(-1)*k.*WPad);
PFRec = pupilPad.*exp(sqrt(-1)*k*WRec);
PFDif=pupilPad.*exp(sqrt(-1)*k*(WRec*0));
PFRes=pupilPad.*exp(sqrt(-1)*k*Residuo);

%PSF
PSFOrig=fft2(PFOrig); %Two-dimensional discrete Fourier Transform.
clear PFOrig;
PSFOrig=fftshift(PSFOrig);
PSFOrig=PSFOrig.*conj(PSFOrig);

PSFRec=fft2(PFRec); %Two-dimensional discrete Fourier Transform.
clear PFRec;
PSFRec=fftshift(PSFRec);
PSFRec=PSFRec.*conj(PSFRec);

PSFDif=fft2(PFDif); %Two-dimensional discrete Fourier Transform.
clear PFDif;
PSFDif=fftshift(PSFDif);
PSFDif=PSFDif.*conj(PSFDif);

PSFRes=fft2(PFRes); %Two-dimensional discrete Fourier Transform.
clear PFRes;
PSFRes=fftshift(PSFRes);
PSFRes=PSFRes.*conj(PSFRes);

                                            if pintar==1
                                                pintaPSF(PSFOrig,PSFRec,PSFRes)
                                                pintarConvolucion(PSFOrig,PSFRec,PSFRes)
                                            end
                                           


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

                                            if pintar==1
                                                pintarMTF(resolucion*2+1,resolucion,LAMBDA,MTFOrig,MTFRec,MTFRes)
                                            end


% *************************************************************************
% *************** QUALITY METRICS**************************************
% *************************************************************************
StRes=max(max(PSFRes))/max(max(PSFDif));
fprintf('Strhel ratio of residual = %2.5f\n',StRes);


StOrig=max(max(PSFOrig))/max(max(PSFDif));
StRec=max(max(PSFRec))/max(max(PSFDif));
fprintf('Difference between Strhel ratio of the incoming phase and the recovered one is %2.5f\n',StOrig-StRec);

RMSOrig=sqrt(sum(c(4:length(c)).^2));
RMSRec=sqrt(sum(ZerRecuperados(4:length(ZerRecuperados)).^2));
RMSResta=sqrt(sum(abs(c(4:length(c))-ZerRecuperados(4:length(ZerRecuperados))).^2));
fprintf('Difference between the RMS of the incoming phase (%2.5f) and the recovered one (%2.5f) is = %2.5f (%2.2f of error in percentage).\n',RMSOrig,RMSRec,RMSResta,100-(RMSRec*100/RMSOrig));
fprintf('The size of the pupil is %2.5f meters\n',TamanoSensor);



                                            figure
                                            plot(4:modos,c(4:modos))
                                            hold on
                                            plot(4:modos,ZerRecuperados(4:modos))
                                            set(gcf,'color','w');
                                            legend('Original','Recovered');
                                            xlabel('Zernike mode')
                                            ylabel('Value in microns')
                                            title('Original Zernike coefficient Vs Recovered')
                                            drawnow();

toc    
