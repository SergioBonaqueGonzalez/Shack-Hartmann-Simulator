function paintPSFs(PSFOrig,PSFRec,PSFRes)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function paints some PSF figures
%}

maximum=max(max(max(PSFOrig)),max(max(PSFRec)));
maximum2=max(maximum,max(max(PSFRes)));
minimum=min(min(min(PSFOrig)),min(min(PSFRec)));
minimum2=min(minimum,min(min(PSFRes)));
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

