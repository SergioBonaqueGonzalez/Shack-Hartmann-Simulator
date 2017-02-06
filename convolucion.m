function [convolucion]=convolucion(PSF)

%Abro el objeto
%[fname,pname] = uigetfile('*.*','Abrir un archivo con el objeto');
myExtract=imread('wooptix.tif');
objeto = myExtract(:,:,1);

%igualo tamaño de matriz
[p1, p2] = size(PSF);
[o1,o2]=size(objeto);

if  rem((p1-o1),2)==0
    objeto2=padarray(objeto,[((p1-o1)/2),((p2-o2)/2)]);
else
    objeto(:,end)=[];
    objeto(end,:)=[];
    [o1,o2]=size(objeto);
    objeto2=padarray(objeto,[((p1-o1)/2),((p2-o2)/2)]);
end


%Hago el producto en el espacio de frecuencias
product =fftshift(fft2(double(objeto2))) .* fftshift(fft2(double(PSF)));

%Vuelvo al espacio imagen
convolution = abs(fftshift(ifft2(double(product))));
clear product;

convolucion = convolution./max(max(convolution));

