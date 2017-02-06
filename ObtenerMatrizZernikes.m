function [VectorZernikes] = ObtenerMatrizZernikes(TamPixel,focalML,Propagacion,lambda,k,MicroLentesMascara,radioMLpxs,CoorBuenas,PupilaML,NFresnel,CentroideReferencia,Escala,modos,ZerModo,Eliminadas,longitud,factor,bits)

Paralelo=parpool('CAFADIS',20);

%Funcion que calcula cada modo de zernike visto a traves del aberrometro
VectorZernikes=zeros(longitud,modos);
PSFmax=zeros(1,modos);
parfor w=1:modos;
    ZMODO=ZerModo{w};
    W=ZMODO.*factor;%paso a metros
    WFSubpupil=W.*MicroLentesMascara;%Frente de onda visto a travez de las microlentes
    Lente=cell(1,length(CoorBuenas));%Prealoco los arrays
    PadLente=cell(1,length(CoorBuenas));%Prealoco los arrays
    
    for i=1:length(CoorBuenas)
        Lente{i}=WFSubpupil(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4));%Separo el frente de onda que "ve" cada microlente. Está en metros
        PadLente{i}=zeros(radioMLpxs*Escala);%Creo una matriz mas grande (escala/2 veces el radio)
        %PadLente{i}(((radioMLpxs*Escala)/2)-radioMLpxs+0.5:end+radioMLpxs-((radioMLpxs*Escala)/2)-0.5,((radioMLpxs*Escala)/2)-radioMLpxs+0.5:end+radioMLpxs-((radioMLpxs*Escala)/2)-0.5)= Lente{i};%Meto la lente en una matriz el doble de grande para eliminar efectos de borde
        PadLente{i}(radioMLpxs*2+1:end-radioMLpxs*2,radioMLpxs*2+1:end-radioMLpxs*2)= Lente{i};%Meto la lente en una matriz el doble de grande para eliminar efectos de borde
    end

    %CALCULO DE LA PSF DE CADA MICROLENTE
    if Propagacion==0 %CASO DE NO PROPAGACION ENTRE MICROLENTES Y CCD
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
            PSF=cell(1,length(PadLente));
            PSFmaxLocal=zeros(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensidad de la imagen 1
                PSFmaxLocal(i)=max(max(PSF{i}));
            end
        elseif NFresnel>= 1
            PSF=cell(1,length(PadLente));
            PSFmaxLocal=zeros(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensidad de la imagen 1
                PSFmaxLocal(i)=max(max(PSF{i}));
            end
        else
            PSF=cell(1,length(PadLente));
            PSFmaxLocal=zeros(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFraunhofer(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensidad de la imagen 1
                PSFmaxLocal(i)=max(max(PSF{i}));
            end
        end
    end
    
    PSFmax(w)=max(max(PSFmaxLocal));
end

PSFmaxTotal=max(max(PSFmax));


parfor w=1:modos;
    ZMODO=ZerModo{w};
    W=ZMODO.*factor;%paso a metros
    WFSubpupil=W.*MicroLentesMascara;%Frente de onda visto a travez de las microlentes
    Lente=cell(1,length(CoorBuenas));%Prealoco los arrays
    PadLente=cell(1,length(CoorBuenas));%Prealoco los arrays
    for i=1:length(CoorBuenas)
        Lente{i}=WFSubpupil(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4));%Separo el frente de onda que "ve" cada microlente. Está en metros
        PadLente{i}=zeros(radioMLpxs*Escala);%Creo una matriz mas grande (escala/2 veces el radio)
        %PadLente{i}(((radioMLpxs*Escala)/2)-radioMLpxs+0.5:end+radioMLpxs-((radioMLpxs*Escala)/2)-0.5,((radioMLpxs*Escala)/2)-radioMLpxs+0.5:end+radioMLpxs-((radioMLpxs*Escala)/2)-0.5)= Lente{i};%Meto la lente en una matriz el doble de grande para eliminar efectos de borde
        PadLente{i}(radioMLpxs*2+1:end-radioMLpxs*2,radioMLpxs*2+1:end-radioMLpxs*2)= Lente{i};%Meto la lente en una matriz el doble de grande para eliminar efectos de borde
    end

    %CALCULO DE LA PSF DE CADA MICROLENTE
    if Propagacion==0 %CASO DE NO PROPAGACION ENTRE MICROLENTES Y CCD
        PSF=cell(1,length(PadLente));%Prealoco los arrays
        for i=1:length(PadLente)
            PF=PupilaML.*exp(sqrt(-1)*k.*PadLente{i});%Funcion pupila
            PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;%PSF
        end
    else %SUPONIENDO PROPAGACION ENTRE MICROLENTE Y CCD
        [~,a]=size(PadLente{1});
        L=a*TamPixel;%Tamaño del lado de cada region considerada
        if NFresnel> 0.5 && NFresnel<1;
            PSF=cell(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensidad de la imagen 1
            end
        elseif NFresnel>= 1
            PSF=cell(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFresnel(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensidad de la imagen 1
            end
        else
            PSF=cell(1,length(PadLente));
            for i=1:length(PadLente)
                S = exp(k*1i.*PadLente{i}); %complex phase screen
                propagada=propagacionFraunhofer(S.*PupilaML,L,lambda,focalML);
                PSF{i} = (abs(propagada).^2);%Intensidad de la imagen 1
            end
        end
    end
    
%     if bits==0
%         for i=1:length(PSF)
%             PSF{i}=PSF{i}/PSFmaxTotal;
%         end
%     else
%         for i=1:length(PSF)
%             PSF{i}=round((PSF{i}/PSFmaxTotal)*((2^bits)-1));
%         end
%     end
    
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
        CoorCentroide(i,1)=CentroideReferencia(i,1)+ LongitudCentroide(i,1);%coordenadas del centroide real
        CoorCentroide(i,2)=CentroideReferencia(i,3)+LongitudCentroide(i,2);
    end
    
    deltas=double(LongitudCentroide*TamPixel); %Estos deltas estan muy discretizados segun el número de pixeles.
    
    alfax= double(atan(deltas(:,1)/focalML));
    alfay=double(atan(deltas(:,2)/focalML));
    
    %pongo cada alfa en su lugar de la matriz inicial
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
    
    %Las junto en un vector. 
    solucionbi=zeros(length(solucionx)*2,1);
    solucionbi(1:length(solucionx),1)=solucionx;
    solucionbi(length(solucionx)+1:length(soluciony)*2,1)=soluciony;
        
    VectorZernikes(:,w)=solucionbi;
end

delete(Paralelo);

