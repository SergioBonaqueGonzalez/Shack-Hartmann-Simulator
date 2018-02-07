function [convolution]=convolution(PSF)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function convolves a PSF with and extended object. 

INPUTS: 
    PSF: the point spread function.

OUTPUTS:
    Convolucion: the convolution of the PSF with the object.
%}
%Open an image
myExtract=imread('wooptix.tif');
objeto = myExtract(:,:,1);

%Resize matrix size
[p1, p2] = size(PSF);
[o1,o2]=size(objeto);
if p1>o1
    if  rem((p1-o1),2)==0
        objeto=padarray(objeto,[((p1-o1)/2),((p2-o2)/2)]);
    else
        objeto(:,end)=[];
        objeto(end,:)=[];
        [o1,o2]=size(objeto);
        objeto=padarray(objeto,[((p1-o1)/2),((p2-o2)/2)]);
    end
elseif p1<o1
    if  rem((o1-p1),2)==0
        objeto(end-((o1-p1)/2)+1:end,:)=[];
        objeto(:,end-((o1-p1)/2)+1:end)=[];
        objeto(:,1:((o1-p1)/2))=[];
        objeto(1:((o1-p1)/2),:)=[];
    else
        objeto(:,end)=[];
        objeto(end,:)=[];
        [o1,o2]=size(objeto);
        objeto(end-((o1-p1)/2)+1:end,:)=[];
        objeto(:,end-((o1-p1)/2)+1:end)=[];
        objeto(:,1:((o1-p1)/2))=[];
        objeto(1:((o1-p1)/2),:)=[];
    end
end

%Make the product in the frequencies space
product =fftshift(fft2(double(objeto))) .* fftshift(fft2(double(PSF)));

%Coming back to the image space.
convolution = abs(fftshift(ifft2(double(product))));
clear product;

convolution = convolution./max(max(convolution));
