function pintarConvolucion(PSFOrig,PSFRec,PSFRes)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This paints some convolution images

 convolucionOrig=convolucion(PSFOrig);
 convolucionRec=convolucion(PSFRec);
 convolucionRes=convolucion(PSFRes);
aux=length(convolucionOrig)/3;
aux2=length(convolucionOrig)-aux;
 
figure;
set(gcf,'color','w');
suptitle('Convolution of PSF with an image') 
subplot(1,3,1);
imshow(convolucionOrig);
title('Original')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])

subplot(1,3,2);
imshow(convolucionRec);
title('Recovered')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])

subplot(1,3,3);
imshow(convolucionRes);
title('Residual')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])
drawnow();
