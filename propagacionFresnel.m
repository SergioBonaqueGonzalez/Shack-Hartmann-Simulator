function[u2]=propagacionFresnel(u1,L,lambda,z,contador)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function calculates the propagation in the Fresnel regime. The transfer function assumes an uniform sampling. 


%INPUTS:
% u1 - Origin field
% L1 - Length of the side of the source 
% lambda - wavelength in meters
% z - propagation distance in meters. 
%contador= dummy variable which indicates if a message should be showed or not

%OUTPUTS:
% L2 - Length of the side of the observed field.  
% u2 - observed field

[M,~]=size(u1); 
dx=L/M; 
crit = abs(lambda*z/L);

if contador==0; %It is only showed the first time
    if dx > crit
        disp('Propagation in Fresnel regime was succesful!');
    else
        disp('Insufficient sampling');
    end
end

fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords 
[FX,FY]=meshgrid(fx,fx); 

H=exp(-1i*pi*lambda*z*(FX.^2+FY.^2)); %trans func.  The exp(jkz) term is ignored. This term doesn’t affect the transverse spatial structure of the observation plane result.
H=fftshift(H); %shift trans func 
U1=fft2(fftshift(u1)); %shift, fft src field 
U2=H.*U1; %multiply 
u2=ifftshift(ifft2(fftshift(U2))); %inv fft, center obs field 
end
