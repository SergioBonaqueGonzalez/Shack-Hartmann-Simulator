function[u2,L2]=propagacionFraunhofer(u1,L1,lambda,z,contador)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function calculates the propagation in the Fraunhofer region. A uniform sampling is supossed. 

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
dx1=L1/M;                 %Sampling intervale
k=2*pi/lambda;            

L2=lambda*z/dx1;          
dx2=lambda*z/L1;          
x2=-L2/2:dx2:L2/2-dx2;    


if contador==0; %It is only showed the first time
    disp('Propagation in Fresnel regime was succesful!');
end

[X2,Y2]=meshgrid(x2,x2); 
c=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X2.^2+Y2.^2)); 
u2=c.*ifftshift(fft2(fftshift(u1)))*dx1^2; 
end 
 
