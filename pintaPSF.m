function pintaPSF(PSFOrig,PSFRec,PSFRes)

maximo=max(max(max(PSFOrig)),max(max(PSFRec)));
maximo2=max(maximo,max(max(PSFRes)));
minimo=min(min(min(PSFOrig)),min(min(PSFRec)));
minimo2=min(minimo,min(min(PSFRes)));
aux=length(PSFOrig)/3;
aux2=length(PSFOrig)-aux;


% figure; %pinta la psf de la pupila original
% suptitle('Point Spread Function todas con la misma escala (arriba) y cada una con su escala (abajo)') 
% subplot(2,3,1)
% imshow(PSFOrig)
% caxis manual
% caxis([minimo2 maximo2]);
% xlim([aux aux2])
% ylim([aux aux2])
% title('PSF original');
% xlabel('\theta_(_x_) (rad)');
% ylabel('\theta_(_y_)(rad)');
% colormap(hot)
% 
% subplot(2,3,2)
% imshow(PSFRec)
% caxis manual
% caxis([minimo2 maximo2]);
% xlim([aux aux2])
% ylim([aux aux2])
% title('PSF recuperada');
% xlabel('\theta_(_x_) (rad)');
% ylabel('\theta_(_y_)(rad)');
% colormap(hot)
% 
% subplot(2,3,3)
% imshow(PSFRes)
% caxis manual
% caxis([minimo2 maximo2]);
% xlim([aux aux2])
% ylim([aux aux2])
% title('PSF residuo');
% xlabel('\theta_(_x_) (rad)');
% ylabel('\theta_(_y_)(rad)');
% colormap(hot)
% 
% subplot(2,3,4)
% imshow(PSFOrig,[])
% title('PSF original');
% xlabel('\theta_(_x_) (rad)');
% ylabel('\theta_(_y_)(rad)');
% xlim([aux aux2])
% ylim([aux aux2])
% colormap(hot)
% 
% subplot(2,3,5)
% imshow(PSFRec,[])
% title('PSF recuperada');
% xlabel('\theta_(_x_) (rad)');
% ylabel('\theta_(_y_)(rad)');
% xlim([aux aux2])
% ylim([aux aux2])
% colormap(hot)
% 
% subplot(2,3,6)
% imshow(PSFRes,[])
% title('PSF residuo');
% xlabel('\theta_(_x_) (rad)');
% ylabel('\theta_(_y_)(rad)');
% xlim([aux aux2])
% ylim([aux aux2])
% colormap(hot)
% 
figure
set(gcf,'color','w');
suptitle('Point Spread Function cada una con su escala') 
subplot(1,3,1)
imshow(PSFOrig,[])
title('PSF original');
xlabel('\theta_(_x_) (rad)');
ylabel('\theta_(_y_)(rad)');
xlim([aux aux2])
ylim([aux aux2])
colormap(hot)

subplot(1,3,2)
imshow(PSFRec,[])
title('PSF recuperada');
xlabel('\theta_(_x_) (rad)');
ylabel('\theta_(_y_)(rad)');
xlim([aux aux2])
ylim([aux aux2])
colormap(hot)

subplot(1,3,3)
imshow(PSFRes,[])
title('PSF residuo');
xlabel('\theta_(_x_) (rad)');
ylabel('\theta_(_y_)(rad)');
xlim([aux aux2])
ylim([aux aux2])
colormap(hot)
drawnow();
zoom

