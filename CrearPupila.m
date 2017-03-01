function [pupil,rho]=CrearPupila(resolucion)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
%This function creates a binary circular pupil inside a square matrix of certain size. 

%INPUTS:
%resolucion= size of the square matrix where the pupil will be inscribed.

%OUTPUTS:
%pupil= a binary circular pupil (inside=1, outside=0) inscribed in the requested matrix size.
%rho= having the pupil a radius of 1, rho is a matrix with the hypotenuse values of each pixel. 

xp=linspace(-1,1,resolucion); 
[X,Y]=meshgrid (xp,xp); 
[rho]=sqrt(X.^2+Y.^2); 


pupil=ones(size(rho));
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;% Defining a circle of radius 1. 
            pupil(i,j)=0;
        end;
    end;
end;
