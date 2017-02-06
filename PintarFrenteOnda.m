function [pupilpintar]=PintarFrenteOnda(rho,W,pintar)

pupilpintar=ones(size(rho));
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;%Este es el valor del radio del circulo unidad donde van a estar definidos los zernikes.
            pupilpintar(i,j)=NaN;
        end;
    end;
end;


if pintar==1
    figure;
    subplot(1,3,1)
    imshow(W.*pupilpintar,[])
    title('Frente de onda original(m)') 
    colorbar
    set(gcf,'color','w');
    drawnow();
end


                                       