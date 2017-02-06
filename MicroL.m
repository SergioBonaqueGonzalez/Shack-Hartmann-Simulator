function [MicroLentesMascara,radioML,CoorBuenas,Eliminadas]=MicroL(NLentes,resolucion,SeparacionML,pupila,pupilpintar,W,pintar)

%Funcion de matlab que dibuja microlentes. Se supone que el detector es
%cuadrado y que hay un espacio 0 entre microlentes

%MicroLentesMascara= mascara con las microlentes
%x11,y11,x22,y22= valores donde se localizan las esquinas superiores izquierdas de
%cada microlente
%AnchoCirculoPequeno= Diametro en pixeles de cada microlente
Fuera = 0; %Valor fuera de los circulo
espaciado=floor((SeparacionML/2)+resolucion/NLentes(1)); %Esto es los pixeles de la CCD asignados a cada microlente
Resto=resolucion-NLentes(1)*espaciado;

if  rem(Resto,2)~=0
    error('Error.La resolucion actual no permite una disposicion simetrica de las microlentes. Pruebe a aumentar o disminuir en uno el numero de microlentes o a cambiar la resolucion.')
end

if Resto<0
    espaciado=espaciado-1;
end

radioML=espaciado/2;
if  rem(espaciado,2)==0
    error('Error.El radio de las microlentes resultantes no permite una configuracion simetrica.Pruebe a aumentar o disminuir en uno el numero de microlentes.')
end


% Crea la pupila de cada microlente
xp=linspace(-1,1,espaciado);

CirculoPequeno = ones(espaciado, espaciado);
[X,Y]=meshgrid (xp,xp); %Meshgrid crea una matriz cuyas filas son copias del vector xp, y cuyas columnas son copias del vector xp
[rho]=sqrt(X.^2+Y.^2); %hipotenusa. Matriz que sustituye cada valor por el de su hipotenusa con respecto a su posicion
[a,b]=size(rho); 
for i=(1:a);
    for j=(1:b);
        if rho(i,j) > 1 ;%Este es el valor del radio del circulo unidad donde van a estar definidos los zernikes.
            CirculoPequeno(i,j)=0;
        end;
    end;
end;


% Voy a definir la esquina superior izquierda donde va alojado el cuadrado que encierra cada circulo
secuenciax=(Resto/2)+1:espaciado:resolucion-espaciado+1; %En el otro eje es igual al trabajar con matrices cuadradas


[p,q] = meshgrid(secuenciax, secuenciax);
pares = [p(:) q(:)];
MuchosCirculos = zeros(resolucion, resolucion);

x11=zeros(1,length(pares));
x22=x11;
y11=x11;
y22=x11;

% Creo los circulo uno a uno
for k = 1 : length(pares)
        % Encuentra la esquina donde se va emplazar cada circulo
        x1 = int16(pares(k,1));
        x11(k)=x1;
        y1 = int16(pares(k,2));
        y11(k)=y1;
        x2 = int16(x1 + espaciado - 1);
        x22(k)=x2;
        y2 = int16(y1 + espaciado - 1);
        y22(k)=y2;
        % Añade el circulo pequeño a la imagen inicial.
        MuchosCirculos(y1:y2, x1:x2) = MuchosCirculos(y1:y2,x1:x2) + CirculoPequeno;
end
% Me aseguro que el valor fuera es cero.
MuchosCirculos(MuchosCirculos == 0) = Fuera;



% Multiplico la imagen de la pupila con la de las microlenes para quedarme solo con las centrales
MicroLentesMascara = MuchosCirculos.*pupila;

%Defino todos los circulos y veo si estan completos. Si no lo estan se
%eliminaran
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
CoorBuenas(1,:)=[]; %La primera linea son todo ceros por la manera que construí estos vectores.

if pintar==1
    subplot(1,3,2)
    imshow(MicroLentesMascara);
    title('Mascara de Microlentes');
    set(gcf,'color','w');
    xlabel('pixeles')
    ylabel('pixeles')
    colorbar

    subplot(1,3,3)
    imshow(pupilpintar.*MicroLentesMascara.*W,[])
    title('Matriz de Microlentes superpuesta al frente de onda')
    set(gcf,'color','w');
    xlabel('pixeles')
    ylabel('pixeles')
    drawnow();
    colorbar
    
end
