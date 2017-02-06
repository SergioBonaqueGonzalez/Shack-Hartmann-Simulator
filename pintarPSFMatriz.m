function pintarPSFMatriz(PSF,MicroLentesMascara,radioMLpxs,CoorBuenas,pupilpintar,Escala)

Sobrante=((Escala*radioMLpxs)-radioMLpxs*2)/2;%Esto son los pixeles que se añadieron de mas a cada lado de las microlentes para calcular la PSF
BackGround=zeros(size(MicroLentesMascara)+2*Sobrante);

%Junto todas las PSF en la misma figura
for i=1:length(PSF)
        BackGround(CoorBuenas(i,1):CoorBuenas(i,2)+2*Sobrante, CoorBuenas(i,3):CoorBuenas(i,4)+2*Sobrante)=BackGround(CoorBuenas(i,1):CoorBuenas(i,2)+2*Sobrante, CoorBuenas(i,3):CoorBuenas(i,4)+2*Sobrante)+PSF{i};

end

figure
subplot(1,2,1)
imshow(BackGround,[])
title('Imagen en la CCD (dobles spots incluidos)') 
set(gcf,'color','w');
xlabel('pixeles')
ylabel('pixeles')



%Ahora me quedo con la parte de la PSF que no sobrepase la microlente. Esto
%es solo para pintar
x=(length(PSF{1})/2)-radioMLpxs;%Regla de tres para encontrar la region central de la PSF
for i=1:length(PSF)
    PSF{i}=PSF{i}(x:end-x-1,x:end-x-1);
end

clear BackGround
BackGround=zeros(size(MicroLentesMascara));

for i=1:length(PSF)
    BackGround(CoorBuenas(i,1):CoorBuenas(i,2), CoorBuenas(i,3):CoorBuenas(i,4))=PSF{i};
end


subplot(1,2,2)
imshow(BackGround.*pupilpintar,[])
title('Imagen en la CCD (respetando region de cada microlente). El punto rojo marca el centroide de referencia') 
xlabel('pixeles')
ylabel('pixeles')
hold on
for i=1:length(PSF)
    rectangle('Position', [CoorBuenas(i,1) CoorBuenas(i,3) radioMLpxs*2 radioMLpxs*2],'LineWidth', 0.1, 'EdgeColor', 'b');
    %pause(0.01)
end
hold on
plot(CoorBuenas(:,1)+radioMLpxs,CoorBuenas(:,3)+radioMLpxs,'o','MarkerSize',1,'MarkerEdgeColor','r')

drawnow();