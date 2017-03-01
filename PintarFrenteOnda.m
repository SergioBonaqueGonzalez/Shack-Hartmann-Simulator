function [pupilpintar]=PintarFrenteOnda(rho,W,pintar)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function paint a wavefront inside a pupil and return a pupil with NaN instead of zeros, useful for cosmetic purposes.

pupilpintar=ones(size(rho));
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;
            pupilpintar(i,j)=NaN;
        end;
    end;
end;


if pintar==1
    figure;
    subplot(1,3,1)
    imshow(W.*pupilpintar,[])
    title('Original wavefront (m)') 
    colorbar
    set(gcf,'color','w');
    drawnow();
end


                                       
