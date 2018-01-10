function [pupil4paint]=PaintWavefront(SH,WF)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function paint a wavefront inside a pupil and return a pupil with NaN instead of zeros, useful for cosmetic purposes.
%}

pupil4paint=ones(size(SH.rho));
[a,b]=size(SH.rho); 
for i=(1:a);
    for j=(1:b);
        if SH.rho(i,j) > 1 ;
            pupil4paint(i,j)=NaN;
        end;
    end;
end;


if SH.paint==1
    figure;
    subplot(1,3,1)
    imshow(WF.*pupil4paint,[])
    title('Original wavefront (m)') 
    colorbar
    set(gcf,'color','w');
    drawnow();
end


                                       
