function pintarWFs(W,Wrec,pupilpintar)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function paints the wavefront.

maximo=max(max(max(W)),max(max(Wrec)));
maximo2=max(maximo,max(max(W-Wrec)));
minimo=min(min(min(W)),min(min(Wrec)));
minimo2=min(minimo,min(min(W-Wrec)));

figure
set(gcf,'color','w');
suptitle('Phases. All with the same scale') 
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
title('ERROR (Original phase minus recovered phase)')
colorbar

subplot(2,2,4);
imshow((W-Wrec).*pupilpintar,[])
title('ERROR (Original phase minus recovered phase) AUTO-SCALED')
colorbar
drawnow();
