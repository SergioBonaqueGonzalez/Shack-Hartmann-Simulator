function mascara=crearmascara(puntos)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function creates a circular mask inside a certain square matrix.

xp=linspace(-1,1,puntos); 
[X,Y]=meshgrid (xp,xp); 
[rho]=sqrt(X.^2+Y.^2); 

[a,b]=size(rho); 
mascara=ones(size(rho));
for i=(1:a)
    for j=(1:b);
        if rho(i,j) > 1 ;%Defines the mask with a radius=1
            mascara(i,j)=0;
        end;
    end;
end;
