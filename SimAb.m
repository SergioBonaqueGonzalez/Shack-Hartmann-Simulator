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

%Todo centroide que esté por fuera de cada microlente está invadiendo la
%microlente contigua. En nuestro simulador la solucion es mejor que en la
%realidad ya que al calcular la PSF y el centroide de la misma
%individualmente y sobre un area mayor no tenemos este problema. En la
%realidad esto es un problema de los sensores de Shack-Hartmann que no
%tendremos con el nuestro. Aun así, vamos a simular una situación ideal
%donde tenemos un algoritmo capaz de identificar que PSF viene de cada
%microlente.

% *************************************************************************
% ***********CALCULO DE LA PENDIENTE***************************************
%**************************************************************************
% *************************************************************************  
%Calculo los deltas en la posicion del centroide en cada microlente en
%metros
deltas=double(LongitudCentroide*TamPixel); %Estos deltas estan muy discretizados segun el número de pixeles.

%Calculo el alfa en cada microlente
alfax= double(atan(deltas(:,1)/focalML));
alfay=double(atan(deltas(:,2)/focalML));

% %pongo cada alfa en su lugar de la matriz inicial
% solucionbipre(1:length(solucionx),1)=solucionx;
% solucionbipre(length(solucionx)+1:length(soluciony)*2,1)=soluciony;



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

      
                                             %Para pintar
                                            if pintar==1
                                                % %Construyo las matrices de pendientes
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
                                                suptitle('Pendientes del frente de onda recuperadas')
                                                subplot(1,2,1)
                                                set(gcf,'color','w');
                                                imshow(solucionx2,[])
                                                title('Eje x') 
                                                colorbar

                                                subplot(1,2,2)
                                                imshow(soluciony2,[])
                                                title('Eje y') 
                                                colorbar
                                                drawnow();

                                                clear solucionx2
                                                clear soluciony2
                                            end

%Las junto en un vector. 
solucionbi=zeros(length(solucionx)*2,1);
solucionbi(1:length(solucionx),1)=solucionx;
solucionbi(length(solucionx)+1:length(soluciony)*2,1)=soluciony;

                                            

% *************************************************************************
% ***********RECUPERACION DE LOS ZERNIKES**********************************
%**************************************************************************
% *************************************************************************  
%El numero mas alto de Zernikes que pueden ser usados en la recuperacion
%está limitado por el numero de microlentes en una direccion. Esto
%significa que si tenemos 30x30 microlentes, no deberiamos exceder el modo 900, en caso contrario el sistema estará sobredeterminado.

%Si quisiera sacar la pendiente del frente de onda original
%[WX,WY] = gradient(W,TamPixel); %Aparece el problema de los bordes y
%habria que quitarlos para verlo bien
% pendiente = sqrt(WX.^2 + WY.^2);
  
%Para la recuperacion de la fase voy a utilizar estimacion lineal sin ligaduras,
%utilizando la variante de estimacion modal. Se basa en la estimación de la aberracion de onda mediante la estimación de los coeficientes modales del desarrollo de la aberración en una base de polinomios ortogonales
%Un buen libro es este: "Wavefront Optics for Vision Correction", Escrito por Guang-ming Dai.
%"https://books.google.es/books?id=aCC-IciqKkYC&pg=PA123&lpg=PA123&dq=matlab+wave+front+derivatives&source=bl&ots=5tSBxQAEbN&sig=ToS1Ns-zQSb6q43NsA02nA513Os&hl=es&sa=X&ved=0ahUKEwiO8fSHkYPRAhUBXRoKHUdQDMwQ6AEINzAD#v=onepage&q=matlab%20wave%20front%20derivatives&f=false"




% ****************CREA EL VECTOR DE ZERNIKES PARA COMPARAR*****************
%Como este proceso tarda mucho, primero se comprueba si ya existe la estructura con las caracteristicas del sensor que se hayan puesto, sino existen si se
%crean

if exist ('VectorZernikes.mat','file') ~=0
    est=load ('VectorZernikes.mat');
    configpre=est.configuracion;
    clear est
    if configpre.resolucion==resolucion && configpre.TamPixel==TamPixel && configpre.focalML == focalML && configpre.Propagacion==Propagacion && configpre.lambda==lambda && configpre.NLentes == NLentes(1) && configpre.modos==modos && configpre.Escala==Escala && configpre.factor==factor && configpre.bits==bits
        fprintf ('Vector de zernikes ya creado para esta configuracion de shack-hartmann.\n')
        VectorZernikes=configpre.VectorZernikes;
    else
        delete('VectorZernikes.mat');
        fprintf ('Se ha modificado la configuracion.Calculando vector de zernikes...\n')
        [VectorZernikes] = ObtenerMatrizZernikes(TamPixel,focalML,Propagacion,lambda,k,MicroLentesMascara,radioMLpxs,CoorBuenas,PupilaML,NFresnel,CentroideReferencia,Escala,modos,ZerModo,Eliminadas,length(solucionbi),factor,bits);
        configuracion = struct('resolucion',resolucion,'TamPixel',TamPixel,'focalML',focalML,'Propagacion',Propagacion,'lambda',lambda,'NLentes',NLentes(1),'modos',modos,'VectorZernikes',VectorZernikes,'Escala',Escala,'factor',factor,'bits',bits);%#ok<NASGU>
        save('VectorZernikes.mat', 'configuracion')
    end
else
    fprintf('Vector de zernikes no disponible.Calculando vector de zernikes...\n')
    [VectorZernikes] = ObtenerMatrizZernikes(TamPixel,focalML,Propagacion,lambda,k,MicroLentesMascara,radioMLpxs,CoorBuenas,PupilaML,NFresnel,CentroideReferencia,Escala,modos,ZerModo,Eliminadas,length(solucionbi),factor,bits);
    configuracion = struct('resolucion',resolucion,'TamPixel',TamPixel,'focalML',focalML,'Propagacion',Propagacion,'lambda',lambda,'NLentes',NLentes(1),'modos',modos,'VectorZernikes',VectorZernikes,'Escala',Escala,'factor',factor,'bits',bits);%#ok<NASGU>                                 
    save('VectorZernikes.mat', 'configuracion')
end

clear configuracion

%Recuperacion por mínimos cuadrados. Es equivalente a construir la matriz
%de recuperacion
ZerRecuperados = lsqr((VectorZernikes(:,1:modos)*1e-6)/factor,solucionbi,1e-10,500);


%frente de onda recuperado
ZSUMARec=zeros(length(ZerModo{1}));
for i=4:modos
    ZMODORec=ZerRecuperados(i)*ZerModo{i};
    ZSUMARec=ZSUMARec+ZMODORec;
end


                                                pintarWFs(W,ZSUMARec.*1e-6,pupilpintar);
                                                
%Voy a poner un borde alrededor de la pupila para hacer la psf
ZTOTALRec=zeros(resolucion*2+1);
ZTOTALPad=ZTOTALRec;
pupilPad=ZTOTALRec;
ZTOTALRec((resolucion/2)+1.5:end-(resolucion/2)-0.5,(resolucion/2)+1.5:end-(resolucion/2)-0.5)=ZSUMARec;
ZTOTALPad((resolucion/2)+1.5:end-(resolucion/2)-0.5,(resolucion/2)+1.5:end-(resolucion/2)-0.5)=ZSUMA;
pupilPad((resolucion/2)+1.5:end-(resolucion/2)-0.5,(resolucion/2)+1.5:end-(resolucion/2)-0.5)=pupil;


WRec=ZTOTALRec.*1e-6;%paso a metros
WPad=ZTOTALPad.*1e-6;%paso a metros
                      

% *************************************************************************
% ***********ANALISIS DEL RESIDUO******************************************
%**************************************************************************
% *************************************************************************                                            
%MATRIZ DE ACOPLAMIENTO. Pagina 69 tesis de Justo Arines
%B es la matriz matriz de deviradas de Zernikes evaluadas en los puntos de medida pero con el numero de Zernikes que tiene el haz incidente
%AA=inv((transpose(VectorZernikes)*VectorZernikes))*transpose(VectorZernikes)*B;
Residuo=pupilPad.*(WPad-WRec);

%Funcion pupila
PFOrig = pupilPad.*exp(sqrt(-1)*k.*WPad);
PFRec = pupilPad.*exp(sqrt(-1)*k*WRec);
PFDif=pupilPad.*exp(sqrt(-1)*k*(WRec*0));
PFRes=pupilPad.*exp(sqrt(-1)*k*Residuo);

%PSF
PSFOrig=fft2(PFOrig); %Two-dimensional discrete Fourier Transform.
clear PFOrig;
PSFOrig=fftshift(PSFOrig);% Esto sería la amplitud compleja de la PSF. Sacado de la ayuda de MATLAB: FFTSHIFT is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
PSFOrig=PSFOrig.*conj(PSFOrig);%Esto ya es el modulo de la PSF

PSFRec=fft2(PFRec); %Two-dimensional discrete Fourier Transform.
clear PFRec;
PSFRec=fftshift(PSFRec);% Esto sería la amplitud compleja de la PSF. Sacado de la ayuda de MATLAB: FFTSHIFT is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
PSFRec=PSFRec.*conj(PSFRec);%Esto ya es el modulo de la PSF

PSFDif=fft2(PFDif); %Two-dimensional discrete Fourier Transform.
clear PFDif;
PSFDif=fftshift(PSFDif);% Esto sería la amplitud compleja de la PSF. Sacado de la ayuda de MATLAB: FFTSHIFT is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
PSFDif=PSFDif.*conj(PSFDif);%Esto ya es el modulo de la PSF

PSFRes=fft2(PFRes); %Two-dimensional discrete Fourier Transform.
clear PFRes;
PSFRes=fftshift(PSFRes);% Esto sería la amplitud compleja de la PSF. Sacado de la ayuda de MATLAB: FFTSHIFT is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
PSFRes=PSFRes.*conj(PSFRes);%Esto ya es el modulo de la PSF

                                            if pintar==1
                                                pintaPSF(PSFOrig,PSFRec,PSFRes)
                                                pintarConvolucion(PSFOrig,PSFRec,PSFRes)
                                            end
                                           


%OTF Y MTF 
OTFOrig=fft2(PSFOrig);%primero calculo la OTF, no hace falta normalizar porque la PSF ya esta normalizada. Esta OTF no esta centrada
MTFOrig = abs(OTFOrig); % MTF es la magnitud (abs) de la OTF
MTFOrig=fftshift(MTFOrig); %centro la MTF
MTFOrig = MTFOrig./max(max(MTFOrig)); % Escalo la MTF para que el pico sea 1

OTFRec=fft2(PSFRec);%primero calculo la OTF, no hace falta normalizar porque la PSF ya esta normalizada. Esta OTF no esta centrada
MTFRec = abs(OTFRec); % MTF es la magnitud (abs) de la OTF
MTFRec=fftshift(MTFRec); %centro la MTF
MTFRec = MTFRec./max(max(MTFRec)); % Escalo la MTF para que el pico sea 1

OTFRes=fft2(PSFRes);%primero calculo la OTF, no hace falta normalizar porque la PSF ya esta normalizada. Esta OTF no esta centrada
MTFRes = abs(OTFRes); % MTF es la magnitud (abs) de la OTF
MTFRes=fftshift(MTFRes); %centro la MTF
MTFRes = MTFRes./max(max(MTFRes)); % Escalo la MTF para que el pico sea 1
                                            
                                            if pintar==1
                                                pintarMTF(resolucion*2+1,resolucion,LAMBDA,MTFOrig,MTFRec,MTFRes)
                                            end


% *************************************************************************
% *************** METRICAS DE CALIDAD**************************************
% *************************************************************************
%Para discernir si la calidad el ajuste es buena obtenemos el residuo (fase
%original - fase recuperada). A este residuo le impondremos una razon de
%strhel superior o igual a 0.8 para considerar que se ha recuperado bien.
%Inferior a 0.8 indicaria que se ha recuperado mal.
StRes=max(max(PSFRes))/max(max(PSFDif));
fprintf('La razon de Strhel del residuo es %2.5f\n',StRes);


%Otra opcion sería obtener la razon de sthrel de la fase original y luego
%la de la fase recuperada y normalizar a la de la fase original
StOrig=max(max(PSFOrig))/max(max(PSFDif));
StRec=max(max(PSFRec))/max(max(PSFDif));
fprintf('La diferencia entre la razon de Strhel del frente de onda original y el recuperado es %2.5f\n',StOrig-StRec);

RMSOrig=sqrt(sum(c(4:length(c)).^2));
RMSRec=sqrt(sum(ZerRecuperados(4:length(ZerRecuperados)).^2));
RMSResta=sqrt(sum(abs(c(4:length(c))-ZerRecuperados(4:length(ZerRecuperados))).^2));
fprintf('La diferencia entre RMS del frente de onda original (%2.5f) y el recuperado (%2.5f) es = %2.5f (%2.2f de error en porcentaje).\n',RMSOrig,RMSRec,RMSResta,100-(RMSRec*100/RMSOrig));
fprintf('El tamaño de la pupila es %2.5f metros\n',TamanoSensor);



                                            figure
                                            plot(4:modos,c(4:modos))
                                            hold on
                                            plot(4:modos,ZerRecuperados(4:modos))
                                            set(gcf,'color','w');
                                            legend('Original','recuperado');
                                            xlabel('Modo de Zernike')
                                            ylabel('Valor en micras')
                                            title('Zernikes originales vs recuperados')
                                            drawnow();

toc    
