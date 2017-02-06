function mascara=crearmascara(puntos)

xp=linspace(-1,1,puntos); %ampliar para ampliar el tamaño de la PSF
[X,Y]=meshgrid (xp,xp); %Meshgrid crea una matriz cuyas filas son copias del vector xp, y cuyas columnas son copias del vector xp
[rho]=sqrt(X.^2+Y.^2); %hipotenusa. Matriz que sustituye cada valor por el de su hipotenusa con respecto a su posicion

[a,b]=size(rho); %con  estas lineas estoy dibujando un circulo
mascara=ones(size(rho));
for i=(1:a)
    for j=(1:b);
        if rho(i,j) > 1 ;%Este es el valor del radio del circulo unidad donde van a estar definidos los zernikes.
            mascara(i,j)=0;
        end;
    end;
end;