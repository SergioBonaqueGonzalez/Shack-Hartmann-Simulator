function pintarWFs(W,Wrec,pupilpintar)

maximo=max(max(max(W)),max(max(Wrec)));
maximo2=max(maximo,max(max(W-Wrec)));
minimo=min(min(min(W)),min(min(Wrec)));
minimo2=min(minimo,min(min(W-Wrec)));

figure
set(gcf,'color','w');
suptitle('Fases todas con la misma escala') 
subplot(2,2,1);
imshow(W.*pupilpintar,[])
caxis manual
caxis([minimo2 maximo2]);
title('Original') 
colorbar

subplot(2,2,2);
imshow(Wrec.*pupilpintar,[])
caxis manual
caxis([minimo2 maximo2]);
title('Recuperada') 
colorbar

subplot(2,2,3);
imshow((W-Wrec).*pupilpintar,[])
caxis manual
caxis([minimo2 maximo2]);
title('ERROR (Fase Original - fase recuperada)')
colorbar

subplot(2,2,4);
imshow((W-Wrec).*pupilpintar,[])
title('ERROR (Fase Original - fase recuperada) AUTO-ESCALADA')
colorbar
drawnow();