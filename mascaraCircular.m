function pupil = mascaraCircular(radio,resolucion)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
Another function that creates a circular mask.
%}

xp=linspace(-1,1,resolucion);
[X,Y]=meshgrid (xp,xp); 
[rho]=sqrt(X.^2+Y.^2); 


pupil=ones(size(rho));
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > radio
            pupil(i,j)=0;
        end;
    end;
end;

