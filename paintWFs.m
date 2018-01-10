function paintWFs(WF,RecWF,pupil4paint)
%{ 
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function paints the wavefront.
%}


maximum=max(max(max(WF)),max(max(RecWF)));
maximum2=max(maximum,max(max(WF-RecWF)));
minimum=min(min(min(WF)),min(min(RecWF)));
minimum2=min(minimum,min(min(WF-RecWF)));

figure
set(gcf,'color','w');
suptitle('Phases. All with the same scale') 
subplot(2,2,1);
imshow(WF.*pupil4paint,[])
caxis manual
caxis([minimum2 maximum2]);
title('Original') 
colorbar

subplot(2,2,2);
imshow(RecWF.*pupil4paint,[])
caxis manual
caxis([minimum2 maximum2]);
title('Recovered') 
colorbar

subplot(2,2,3);
imshow((WF-RecWF).*pupil4paint,[])
caxis manual
caxis([minimum2 maximum2]);
title('ERROR (Original-recovered)')
colorbar

subplot(2,2,4);
imshow((WF-RecWF).*pupil4paint,[])
title('ERROR (Original-recovered) AUTO-SCALED')
colorbar
drawnow();
