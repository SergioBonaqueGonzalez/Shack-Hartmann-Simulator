function [pupil,rho]=CrearPupila(resolucion)

% *************************************************************************
% ************************CREACION PUPILA *********************************
% *************************************************************************
% *************************************************************************

%Creacion pupila circular inicial para poder definir zernikes dentro
xp=linspace(-1,1,resolucion); 
[X,Y]=meshgrid (xp,xp); %Meshgrid crea una matriz cuyas filas son copias del vector xp, y cuyas columnas son copias del vector xp
[rho]=sqrt(X.^2+Y.^2); %hipotenusa. Matriz que sustituye cada valor por el de su hipotenusa con respecto a su posicion


pupil=ones(size(rho));
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;%Este es el valor del radio del circulo unidad donde van a estar definidos los zernikes.
            pupil(i,j)=0;
        end;
    end;
end;

%                                         imagesc(Coordenadas,Coordenadas,pupil);%
%                                         title('Pupila Inicial') 
%                                         set(gcf,'color','w');
%                                         xlabel('metros')
%                                         ylabel('metros')

