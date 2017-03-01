function pintarPSFMatriz(PSF,MicroLentesMascara,radioMLpxs,CoorBuenas,pupilpintar,Escala)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function paints the PSF of each microlense.



Sobrante=((Escala*radioMLpxs)-radioMLpxs*2)/2;%These are the virtually added pixels to each microlenses in order to calculate the PSF.
BackGround=zeros(size(MicroLentesMascara)+2*Sobrante);

%Joinning all PSF in the same figure
for i=1:length(PSF)
        BackGround(CoorBuenas(i,1):CoorBuenas(i,2)+2*Sobrante, CoorBuenas(i,3):CoorBuenas(i,4)+2*Sobrante)=BackGround(CoorBuenas(i,1):CoorBuenas(i,2)+2*Sobrante, CoorBuenas(i,3):CoorBuenas(i,4)+2*Sobrante)+PSF{i};
end

figure
subplot(1,2,1)
imshow(BackGround,[])
title('Image in the CCD (doubles spots are included)') 
set(gcf,'color','w');
xlabel('pixels')
ylabel('pixels')



%Now, only portion of PSF which does not exceed the area of each microlense is taken into account. This is only for painting purposes
x=(length(PSF{1})/2)-radioMLpxs;% Find the central area of each microlense.
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
title('Image in the CCD (the area of each microlenses is respected). The red point is the reference centroid') 
xlabel('pixels')
ylabel('pixels')
hold on
for i=1:length(PSF)
    rectangle('Position', [CoorBuenas(i,1) CoorBuenas(i,3) radioMLpxs*2 radioMLpxs*2],'LineWidth', 0.1, 'EdgeColor', 'b');
    %pause(0.01)
end
hold on
plot(CoorBuenas(:,1)+radioMLpxs,CoorBuenas(:,3)+radioMLpxs,'o','MarkerSize',1,'MarkerEdgeColor','r')

drawnow();
