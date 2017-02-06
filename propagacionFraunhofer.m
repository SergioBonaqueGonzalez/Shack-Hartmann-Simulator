function[u2,L2]=propagacionFraunhofer(u1,L1,lambda,z,contador)
% propagacion en campo lejano en el regimen de Fraunhofer. Asume un muestreo uniforme. 
% u1 - Campo fuente/origen
% L1 - Longitud del lado del campo fuente
% lambda - longitud de onda en metros
% z - distancia de propagacion en metros
% L2 - Longitud del lado del campo observado 
% u2 - campo observado

[M,~]=size(u1);           %Obtiene las dimensiones del campo fuente
dx1=L1/M;                 %Intervalo de muestreo
k=2*pi/lambda;            

L2=lambda*z/dx1;          % L2 - Longitud del lado del campo observado 
dx2=lambda*z/L1;          %intervalo del muestreo en el campo observado
x2=-L2/2:dx2:L2/2-dx2;    %coordenadas en el campo observado


if contador==0; %Para que solo lo muestre la primera vez;
    disp('Propagacion en regimen de Fresnel OK!');
end

[X2,Y2]=meshgrid(x2,x2); 
c=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X2.^2+Y2.^2)); 
u2=c.*ifftshift(fft2(fftshift(u1)))*dx1^2; 
end 
 