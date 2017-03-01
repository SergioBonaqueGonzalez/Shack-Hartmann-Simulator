function [convolucion]=convolucion(PSF)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function convolves a PSF with and extended object. 

%INPUTS: 
%PSF: the point spread function.

%OUTPUTS:
%Convolucion: the convolution of the PSF with the object.

%Open an image
myExtract=imread('wooptix.tif');
objeto = myExtract(:,:,1);

%Resize matrix size
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


%Make the product in the frequencies space
product =fftshift(fft2(double(objeto2))) .* fftshift(fft2(double(PSF)));

%Coming back to the image space.
convolution = abs(fftshift(ifft2(double(product))));
clear product;

convolucion = convolution./max(max(convolution));

