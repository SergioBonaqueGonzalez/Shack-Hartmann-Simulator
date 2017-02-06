function[u2]=propagacionFresnel(u1,L,lambda,z,contador)
% propagacion de Fresnel - La funcion de transferencia asume x e y de la
% misma longitud y muestreadas uniformemente.
% u1 - Plano de la fuente
% L - Cuanto mide el lado del plano de observacion y el plano de la fuente
% lambda - longitud de onda en metros
% z - distancia de propagacion
% u2 - plano de observacion
[M,~]=size(u1); %Encuentra el tamaño de la imagen de entrada
dx=L/M; %Intervalo de muestreo
crit = abs(lambda*z/L);
%disp(['dx = ', num2str(dx), ', lambda*z/L = ',num2str(crit)]);
if contador==0; %Para que solo lo muestre la primera vez;
    if dx > crit
        disp('Propagacion en regimen de Fresnel OK!');
    else
        disp('Muestreo insuficiente');
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