function [MicroLentesMascara,radioML,CoorBuenas,Eliminadas]=MicroL(NLentes,resolucion,SeparacionML,pupila,pupilpintar,W,pintar)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function calculates a microlenses array. It is supposed a square detector and no separation between microlenses.

%INPUTS:
%NLentes= number of microlenses in a row.
%Resolucion = number o pixels of the detector
%SeparacionML =  Separation between microlenses. IN THIS MOMENT THIS FEATURE IT IS NOT IMPLEMENTED
%pupila = pupil of the whole system.
%pupilpintar = a dummy pupil with NaN instead of zeros. It is created automatically in a previous script. It is useful for visualization purposes.
%W = Matrix which contains the incoming wavefront.
%pintar = flag that indicates if figures should be painted or not.

%OUTPUTS
%MicroLentesMascara= binary mask which contains the microlenses mask.
%radioML= Radius of each microlense in pixels.
%CoorBuenas= Coordinates of the upper left corner of each valid microlense. A valid microlense is the one which is completely inside the whole pupil of the system.
%Eliminadas = Coordinates of the upper left corner of each invalid microlense. For example, discarded microlenses because they are incomplete.


Fuera = 0; %Assigned value for pixels outside microlenses 
espaciado=floor((SeparacionML/2)+resolucion/NLentes(1)); %Assigned pixels for each microlense.
Resto=resolucion-NLentes(1)*espaciado;

%For testing purposes, i needed a symmetrical disposition of the microlenses array with respect the incoming wavefront:
if  rem(Resto,2)~=0
    error('Error.Actual resolution does not allow a symmetric disposition of microlenses array. Try to increase or diminish the number of microleneses in a row.')
end

if Resto<0
    espaciado=espaciado-1;
end

radioML=espaciado/2;
if  rem(espaciado,2)==0
    error('Error.The radius of microlenses does not allow a symmetric disposition of microlenses array. Try to increase or diminish the number of microleneses in a row.')
end


% Here, the pupil of each microlense will be created.
xp=linspace(-1,1,espaciado);
CirculoPequeno = ones(espaciado, espaciado);
[X,Y]=meshgrid (xp,xp); 
[rho]=sqrt(X.^2+Y.^2); 
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;
            CirculoPequeno(i,j)=0;
        end;
    end;
end;


% the upper left corner of the square where each microlense is inscribed is defined .
secuenciax=(Resto/2)+1:espaciado:resolucion-espaciado+1; 
[p,q] = meshgrid(secuenciax, secuenciax);
pares = [p(:) q(:)];
MuchosCirculos = zeros(resolucion, resolucion);
x11=zeros(1,length(pares));
x22=x11;
y11=x11;
y22=x11;


% Creating each circle one by one:
for k = 1 : length(pares)
        % find upper left corner:
        x1 = int16(pares(k,1));
        x11(k)=x1;
        y1 = int16(pares(k,2));
        y11(k)=y1;
        x2 = int16(x1 + espaciado - 1);
        x22(k)=x2;
        y2 = int16(y1 + espaciado - 1);
        y22(k)=y2;
        % Adding the tiny pupil
        MuchosCirculos(y1:y2, x1:x2) = MuchosCirculos(y1:y2,x1:x2) + CirculoPequeno;
end
MuchosCirculos(MuchosCirculos == 0) = Fuera;% make sure the value outside pupils is zero.

% Multiplication of the main mask and the microlenses mask. 
MicroLentesMascara = MuchosCirculos.*pupila;

% Now, only those microlenses which are completely inside the main pupil are selected: 
CoorBuenas=zeros(1,4);
Eliminadas=zeros(length(y11),1);
for i=1:length(x11)
    if sum(sum(MicroLentesMascara(y11(i):y22(i), x11(i):x22(i))))<sum(sum(CirculoPequeno))
        MicroLentesMascara(y11(i):y22(i), x11(i):x22(i))=0;
    else 
        CoorBuenas(end+1,1:4)=[y11(i),y22(i),x11(i),x22(i)];
        Eliminadas(i)=1;
    end
end
CoorBuenas(1,:)=[]; %In the way these vectors are constructed, the first line is always zero.


%Painting...
if pintar==1
    subplot(1,3,2)
    imshow(MicroLentesMascara);
    title('Microlenses mask');
    set(gcf,'color','w');
    xlabel('pixels')
    ylabel('pixels')
    colorbar

    subplot(1,3,3)
    imshow(pupilpintar.*MicroLentesMascara.*W,[])
    title('Microlenses mask and incoming wavefront')
    set(gcf,'color','w');
    xlabel('pixels')
    ylabel('pixels')
    drawnow();
    colorbar
    
end
