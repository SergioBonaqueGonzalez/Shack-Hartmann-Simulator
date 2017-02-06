function pintarConvolucion(PSFOrig,PSFRec,PSFRes)

 convolucionOrig=convolucion(PSFOrig);
 convolucionRec=convolucion(PSFRec);
 convolucionRes=convolucion(PSFRes);
aux=length(convolucionOrig)/3;
aux2=length(convolucionOrig)-aux;
 
figure;
set(gcf,'color','w');
suptitle('Convolucion PSF con una imagen') 
subplot(1,3,1);
imshow(convolucionOrig);
title('Original')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])

subplot(1,3,2);
imshow(convolucionRec);
title('Recuperada')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])

subplot(1,3,3);
imshow(convolucionRes);
title('Residuo')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])
drawnow();