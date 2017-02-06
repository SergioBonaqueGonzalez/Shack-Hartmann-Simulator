function [ZerRecuperados]=SimAb(modos,c,resolucion,TamPixel,LAMBDA,NLentes,focalML,Propagacion,factor,bits,pintar)

%Programa hecho por Sergio Bonaque-Gonzalez
%   sergio.bonaque@um.es


%para simular un sensor de Shack-Hartmann
%de dos maneras: uno suponiendo propagación del frente de onda entre cada
%microlente y la CCD, por lo que la determinación del centroide es más
%compleja y se acumulan más errores. La segunda manera es sin suponer
%propagación, por lo que al final se trata de hacer un pseudo-binning (microesferas) de la
%imagen, calcular el centroide de cada subapertura, extrapolar las
%pendientes y obtener el frente de onda a partir de las mismas.

%Fecha de inicio: 30/11/2016 

%***modos= modos de zernike que se intentaran recuperar (ej. modos=10)
%***c = vector que contiene los coeficientes de zernike originales en la
%notacion Noll. (ej. c=rand(10,1)
%***resolucion=resolucion de la CCD (se suponen CCD cuadrados e impares) (ej. 1024)
%***TamPixel=Tamaño del pixel en la CCD (ej. 1.471e-6)
%***LAMBDA= longitud de onda en micras (ej. 0.780)
%***NLentes= numero de lentes en cada direccion (se suponen matrices de
%***microlentes cuadradas (ej.41)
%***focalML=Focal de las microlentes en metros
%***Propagacion= marcador que indica si se va a usar propagacion entre las
%microlentes y la CCD: =0 no se considera que haya propagacion entre las lentes y el ccd; =1 SI hay propagacion 
%***factor=el valor con el que se caracterizan los zernikes. ej=1e-8
%***bits=numero de bits que dispone la CCD: si se deja en 0, se realizan los
%calculos con precision doble de MATLAB
%***pintar=Variable que indica si se van a pintar figuras o no 1=SI, 2=NO

%%Ejemplo de uso:
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
c=c(:);
warning('off', 'Images:initSize:adjustingMag'); %Esta linea es para que no muestre el molesto aviso de que se está cambiando el tamaño de la imagen al representar


if  rem(resolucion,2)==0
    error('Error.La resolucion debe ser impar para tener un pixel central y definir los coeficientes de zernike simetricamente.')
end

espaciado=floor(resolucion/NLentes(1)); %Esto es los pixeles de la CCD asignados a cada microlente
Resto=resolucion-NLentes(1)*espaciado;

if  rem(Resto,2)~=0
    error('Error.La configuracion actual no permite una disposicion simetrica de las microlentes. Pruebe a aumentar o disminuir en uno el numero de microlentes o a cambiar la resolucion.')
end



lambda = LAMBDA*1e-6;%Longitud de onda en metros
k = 2*pi/lambda;
TamanoSensor=TamPixel*resolucion;
SeparacionML=0; %Separacion entre microlentes. En este script se supone que es cero. Queda para el futuro implementar una separacion.

Radio=TamanoSensor/2;%Radio de la pupila circunscrita en el sensor.
[pupil,rho]=CrearPupila(resolucion);


% *************************************************************************
% ************************CREACION MATRIZ DE ZERNIKES**********************
% *************************************************************************
% *************************************************************************
%Busca si se ha creado previamente la matriz de zernikes unidad con la
%configuracion actual. 

if exist ('MatrizZernikes.mat','file') ~=0 && exist ('config.mat','file') ~=0
    est=load ('config.mat');
    config=est.config;
    clear est
    if config.resolucion==resolucion && config.modos==modos 
        fprintf ('Matriz de zernikes ya creada.\n')
        ZerModo=load ('MatrizZernikes.mat');
        ZerModo=ZerModo.ZerModo;
    else
        delete('MatrizZernikes.mat');
        delete('config.mat');
        fprintf ('Se ha modificado la configuracion.Calculando nueva matriz de zernikes...\n')
        ZerModo=cell(modos,1);
        for i=1:modos
            ZerModo{i}=zernike(i,resolucion);
        end
        config=struct('resolucion',resolucion,'modos',modos); %#ok<NASGU>
        save('MatrizZernikes.mat', 'ZerModo')
        save('config.mat', 'config')
    end
else
    fprintf('Matriz de zernikes no disponible. Calculando...\n')
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
% *********** INTRODUCIR ZERNIKES EN MICRAS (NOTACION NOLL)****************
% ********** DEL FRENTE DE ONDA SIMULADO A RECUPERAR***********************
%Zernike polyomials and atmospheric turbulence J Op Soc Am. Vol 66, No 3 , March 1976**********************
% *************************************************************************
ZSUMA=zeros(length(ZerModo{1}));

%Voy a suponer que los 3 primeros coeficientes de zernike son cero
for i=2:modos
    ZMODO=c(i)*ZerModo{i};
    ZSUMA=ZSUMA+ZMODO;
end
%Meto los zernikes creados en un cuadrado mas grande como se definió al principio
W=ZSUMA.*1e-6;%paso a metros


                                              

                                            %Para pintar
                                            [pupilpintar]=PintarFrenteOnda(rho,W,pintar);
                                            
                                            
% *************************************************************************
% ***************DEFINICION MATRIZ MICROLENTES*****************************
%**************************************************************************
% *************************************************************************          
[MicroLentesMascara,radioMLpxs,CoorBuenas,Eliminadas]=MicroL(NLentes,resolucion,SeparacionML,pupil,pupilpintar,W,pintar);


% *************************************************************************
% ***********BONDAD DEL MUESTREO EN CADA MICROLENTE************************
%**************************************************************************
% *************************************************************************     
RadioMLm=TamPixel*radioMLpxs;%Radio en mm de las microlentes
FML=focalML/(RadioMLm*2);%Numero F de las microlentes
Deltax=FML*lambda; %Ecuaciones sacadas del libro de Fourier para saber si el muestreo de cada microlente es suficiente
MuestreoNecesario=RadioMLm*2/Deltax;

if MuestreoNecesario<=radioMLpxs*2
    fprintf ('Muestreo de microlentes adecuado.\n')
else
    fprintf ('ERROR: display de microlentes inadecuado.\n')
end

Lente=cell(1,length(CoorBuenas));%Prealoco los arrays
PadLente=cell(1,length(CoorBuenas));%Prealoco los arrays


% *************************************************************************
% ***********CALCULO DE LA IMAGEN DE CADA MICROLENTES**********************
%**************************************************************************
% *************************************************************************     
WFSubpupil=W.*MicroLentesMascara;%Frente de onda visto a traves de las microlentes
%Voy a separar cada microlente para tratarla por separado
Escala=6; %Este valor simboliza el espacio que le dejo a los lados para calcular la PSF. El cuadrado donde va alojada la PSF es Escala veces mas grande. El valor optimo para que no influya en la recuperacion del centroide es el doble.
for i=1:length(CoorBuenas)
        Lente{i}=WFSubpupil(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4));%Separo el frente de onda que "ve" cada microlente. Está en metros
        PadLente{i}=zeros(radioMLpxs*Escala);%Creo una matriz mas grande (escala/2 veces el radio)
        %PadLente{i}(((radioMLpxs*Escala)/2)-radioMLpxs+0.5:end+radioMLpxs-((radioMLpxs*Escala)/2)-0.5,((radioMLpxs*Escala)/2)-radioMLpxs+0.5:end+radioMLpxs-((radioMLpxs*Escala)/2)-0.5)= Lente{i};%Meto la lente en una matriz el doble de grande para eliminar efectos de borde
        PadLente{i}(radioMLpxs*2+1:end-radioMLpxs*2,radioMLpxs*2+1:end-radioMLpxs*2)= Lente{i};%Meto la lente en una matriz el doble de grande para eliminar efectos de borde
end

%Me defino la pupila de cada microlente en binario. La defino fuera de la
%funcion MicroL porque esta pupila ocupa mas pixeles para muestrear bien la
%PSF
PupilaML=PadLente{1};
[a,b]=size(PupilaML);
for i=1:a
    for j=1:b
        if PupilaML(i,j)~=0
           PupilaML(i,j)=1;
        end
    end
end

%CALCULO DE LA PSF DE CADA MICROLENTE
NFresnel=((TamanoSensor/2)^2)/(lambda*focalML); %Si este numero es mucho menor que 1, se usa Fraunhofer. Si no, se puede usar fresnel. Fresnel expression describes diffraction under the paraxial assumption, where only rays that make a small angle (< ~0.1 rad) relative to the optical axis are considered. 
if Propagacion==0 %CASO DE NO PROPAGACION ENTRE MICROLENTES Y CCD
    fprintf('No se considera propagacion entre microlentes y CCD.\n')
    PSF=cell(1,length(PadLente));%Prealoco los arrays
    PSFmaxLocal=zeros(1,length(PadLente));
    for i=1:length(PadLente)
        PF=PupilaML.*exp(sqrt(-1)*k.*PadLente{i});%Funcion pupila
        PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;%PSF
        PSFmaxLocal(i)=max(max(PSF{i}));
    end
else %SUPONIENDO PROPAGACION ENTRE MICROLENTE Y CCD
    [~,a]=size(PadLente{1});
    L=a*TamPixel;%Tamaño del lado de cada region considerada
    if NFresnel> 0.5 && NFresnel<1;
        fprintf('Numero de Fresnel =%2.0f.\n',NFresnel);
        fprintf('Cuidado, estamos en propagacion en el límite del regimen de Fresnel.\n')
        PSF=cell(1,length(PadLente));
        PSFmaxLocal=zeros(1,length(PadLente));
        for j=1:length(PadLente)
            S = exp(k*1i.*PadLente{j}); %complex phase screen
            propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
            PSF{j} = (abs(propagada).^2);%Intensidad de la imagen 1
            PSFmaxLocal(j)=max(max(PSF{j}));
        end
    elseif NFresnel>= 1
        fprintf('Numero de Fresnel =%2.0f.\n',NFresnel);
        fprintf ('Propagacion en regimen de Fresnel.\n')
        PSF=cell(1,length(PadLente));
        PSFmaxLocal=zeros(1,length(PadLente));
        for j=1:length(PadLente)
            S = exp(k*1i.*PadLente{j}); %complex phase screen
            propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
            PSF{j} = (abs(propagada).^2);%Intensidad de la imagen 1
            PSFmaxLocal(j)=max(max(PSF{j}));
        end
    else
        fprintf('Numero de Fresnel =%2.0f.\n',NFresnel);
        fprintf ('Propagacion en regimen de Fraunhofer.\n')
        PSF=cell(1,length(PadLente));
        PSFmaxLocal=zeros(1,length(PadLente));
        for j=1:length(PadLente)
            S = exp(k*1i.*PadLente{j}); %complex phase screen
            propagada=propagacionFraunhofer(S.*PupilaML,L,lambda,focalML);
            PSF{j} = (abs(propagada).^2);%Intensidad de la imagen 1
            PSFmaxLocal(j)=max(max(PSF{j}));
        end
    end
end

%Lo que he hecho para simular una CCD real es normalizar con respecto al
%maximo total. Para simular la profundidad de bits, discretizo con respecto
%al numero de bits
if bits==0
    fprintf('Calculos realizados con precision "double" de MATLAB. \n');
    for j=1:length(PSF)
        PSF{j}=PSF{j}/max(max(PSF{j}));
    end
else
    
    PSFmax=max(max(PSFmaxLocal));
    fprintf('Calculos realizados para %2.0f bits. \n',bits);
    for i=1:length(PSF)
        PSF{i}=round((PSF{i}/PSFmax)*((2^bits)-1));
    end
end


                                            %Para pintar
                                            if pintar==1
                                                pintarPSFMatriz(PSF,MicroLentesMascara,radioMLpxs,CoorBuenas,pupilpintar,Escala)
                                                drawnow();
                                            end
            
% *************************************************************************
% ***********CALCULO DEL CENTROIDE*****************************************
%**************************************************************************
% *************************************************************************   
%Hallar el centroide en un sensor de Shack-Hartmann es un tema complicado.
%El ruido y el numero de fotones es importante y puede influir
%significativamente en el algoritmo que se utilice para encontrar el
%centroide. Cuando estemos con un prototipo real habra que tener en cuenta
%todos estos factores y a partir de la medida de lentes calibradas elegir
%el mejor método de detección del algoritmo. En el caso de hacer
%simulaciones, vamos a aproximar que el detector no genera ruido y que no
%se producen artefactos espureos, por lo que podemos calcular sin temor a
%equivocarnos el centroide directamente con el toolbox de matlab.Para ver
%algunas consecuencias de este problema leer: http://www.ctio.noao.edu/soar/sites/default/files/SAM/archive/5490-123.pdf
%Y algunos métodos para mejorar el calculo en un caso real en: "Shack-Hartmann wavefront sensor image analysis: a comparison of centroiding methods and image-processing techniques"

%La PSF ya esta normalizada del paso anterior
centroideRealX=zeros(1,length(PSF));
centroideRealY=zeros(1,length(PSF));

%Llamo a la funcion CalCentroides que calcula la posicion del centroide con
%precision subpixel.
%En la realidad, existen interferencias entre microlentes que crean un
%segundo o mas spots secundarios que afectan al calculo, siendo el error
%introducido debido a esto =(theta(t)-theta'(t))/theta(t); siendo tetha el
%angulo real de tilt y theta' el calculado ("Analysis on ShackHartmann
%wave-front sensor with Fourier optics"). No obstante, en nuestro caso he
%simulado un sistema perfecto donde cada microlente forma su imagen por
%separado y no tendre en cuenta esto.
display=0; %1=se pintan los centroides, 0=no se pintan
if display==1
    figure
end
for i=1:length(PSF)
    [centroideRealX(i),centroideRealY(i)]=CalCentroides(PSF{i},1,display);
end


%Hallo las coordenadas de los centroides y de los centroides de referencia
CentroideReferencia=CoorBuenas+radioMLpxs;%coordenadas del centroide de referencia
errores=0;
LongitudCentroide=zeros(length(CoorBuenas),2);
CoorCentroide=zeros(length(CoorBuenas),2);
for i=1:length(CentroideReferencia)
    LongitudCentroide(i,1)=centroideRealX(i)-(Escala*radioMLpxs/2);
    LongitudCentroide(i,2)=centroideRealY(i)-(Escala*radioMLpxs/2);
    CoorCentroide(i,1)=CentroideReferencia(i,1)+ LongitudCentroide(i,1);%coordenadas del centroide real
    CoorCentroide(i,2)=CentroideReferencia(i,3)+LongitudCentroide(i,2);
    %Voy a poner un contador para saber cuando el centroide estaría en el
    %area de otra microlente
    if LongitudCentroide(i,1)>radioMLpxs
        errores=errores+1;
    else
        if LongitudCentroide(i,2)>radioMLpxs
            errores=errores+1;
        end
    end
end


                                            %Para pintar
                                            if pintar==1
                                                figure
                                                suptitle('Centroide de referencia vs centroide encontrado')
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
                                                legend('Centroide de referencia','Centroide encontrado')
                                                title('Posicion de los centroides') 
                                                xlabel({'pixeles';['Numero de centroides que se saldrian del area de cada microlente (crosstalk): ' num2str(errores)]});
                                                ylabel('pixeles')
                                                set(gca,'ydir','reverse')
                                                xlim([0 resolucion])
                                                ylim([0 resolucion])
                                                
                                                subplot(1,2,2)
                                                quiver(CentroideReferencia(:,1),CentroideReferencia(:,3),LongitudCentroide(:,1),LongitudCentroide(:,2),0);
                                                title('Desplazamiento del centroide real') 
                                                xlabel('pixeles')
                                                ylabel('pixeles')
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
