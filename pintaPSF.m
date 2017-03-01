function pintaPSF(PSFOrig,PSFRec,PSFRes)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function paints some PSF figures

maximo=max(max(max(PSFOrig)),max(max(PSFRec)));
maximo2=max(maximo,max(max(PSFRes)));
minimo=min(min(min(PSFOrig)),min(min(PSFRec)));
minimo2=min(minimo,min(min(PSFRes)));
aux=length(PSFOrig)/3;
aux2=length(PSFOrig)-aux;


figure
set(gcf,'color','w');
suptitle('Point Spread Function auto-scaled') 
subplot(1,3,1)
imshow(PSFOrig,[])
title('Original PSF');
xlabel('\theta_(_x_) (rad)');
ylabel('\theta_(_y_)(rad)');
xlim([aux aux2])
ylim([aux aux2])
colormap(hot)

subplot(1,3,2)
imshow(PSFRec,[])
title('Recovered PSF');
xlabel('\theta_(_x_) (rad)');
ylabel('\theta_(_y_)(rad)');
xlim([aux aux2])
ylim([aux aux2])
colormap(hot)

subplot(1,3,3)
imshow(PSFRes,[])
title('Residual PSF');
xlabel('\theta_(_x_) (rad)');
ylabel('\theta_(_y_)(rad)');
xlim([aux aux2])
ylim([aux aux2])
colormap(hot)
drawnow();
zoom

