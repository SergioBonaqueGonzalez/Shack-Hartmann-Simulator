function PaintConvolution(PSFOrig,PSFRec,PSFRes)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This paints some convolution images
%}

convolutionOrig=convolution(PSFOrig);
convolutionRec=convolution(PSFRec);
convolutionRes=convolution(PSFRes);
aux=length(convolutionOrig)/3;
aux2=length(convolutionOrig)-aux;

figure;
set(gcf,'color','w');
suptitle('Convolution of PSF with an image')
subplot(1,3,1);
imshow(convolutionOrig);
title('Original')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])

subplot(1,3,2);
imshow(convolutionRec);
title('Recovered')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])

subplot(1,3,3);
imshow(convolutionRes);
title('Residual')
colormap(gray);
xlim([aux aux2])
ylim([aux aux2])
drawnow();
