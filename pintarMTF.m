function pintarMTF(M,L,LAMBDA,MTFOrig,MTFRec,MTFRes)

%%PINTADO Y ESCALADO DE LA MTF
porcion = M/(L/2);%pixeles que corresponden a la funcion pupila
cutofffull = (0.95*2)/(LAMBDA/1000)/57.3;
cutoffaper = (0.95*2)/(LAMBDA/1000)/57.3;
halfsize = ceil(0.5*(cutoffaper/cutofffull)*porcion);
cutoff=cutofffull/2;

% definicion ejes en ciclos/grado
n = size(MTFOrig);
axisMTF = -cutoff:2*cutoff/n(1):(cutoff-(2*cutoff/n(1)));

%MTF radial
[sizMTFOrig, ~] = size(MTFOrig);
[~, peakyOrig]=max(max(MTFOrig));
[~, peakxOrig]=max(max(transpose(MTFOrig)));

[sizMTFRec, ~] = size(MTFRec);
[~, peakyRec]=max(max(MTFRec));
[~, peakxRec]=max(max(transpose(MTFRec)));

[sizMTFRes, ~] = size(MTFRes);
[~, peakyRes]=max(max(MTFRes));
[~, peakxRes]=max(max(transpose(MTFRes)));


maxradiusOrig=sqrt(2*(sizMTFOrig/2)^2);
maxradiusRec=sqrt(2*(sizMTFRec/2)^2);
maxradiusRes=sqrt(2*(sizMTFRes/2)^2);

numbins=100;%sizMTF/4; 
radialMTFarrayOrig = zeros(numbins,2);
radialMTFarrayRec = zeros(numbins,2);
radialMTFarrayRes = zeros(numbins,2);

for i=1:sizMTFOrig
   for j=1:sizMTFOrig
      xposOrig=i-peakxOrig;
      yposOrig=j-peakyOrig;
      radiusOrig=ceil(numbins*sqrt(xposOrig^2+yposOrig^2)/maxradiusOrig);
      if (radiusOrig==0)
      else
         radialMTFarrayOrig(radiusOrig,1)=radialMTFarrayOrig(radiusOrig,1)+MTFOrig(i,j);
         radialMTFarrayOrig(radiusOrig,2)=radialMTFarrayOrig(radiusOrig,2)+1;
      end
   end
end


radialMTFOrig=radialMTFarrayOrig(:,1)./radialMTFarrayOrig(:,2);
radialMTFOrig=cat(1,1,radialMTFOrig);
axisradMTFOrig=0:sqrt(2)*max(abs(axisMTF))/numbins:sqrt(2)*max(abs(axisMTF));

for i=1:sizMTFRec
   for j=1:sizMTFRec
      xposRec=i-peakxRec;
      yposRec=j-peakyRec;
      radiusRec=ceil(numbins*sqrt(xposRec^2+yposRec^2)/maxradiusRec);
      if (radiusRec==0)
      else
         radialMTFarrayRec(radiusRec,1)=radialMTFarrayRec(radiusRec,1)+MTFRec(i,j);
         radialMTFarrayRec(radiusRec,2)=radialMTFarrayRec(radiusRec,2)+1;
      end
   end
end


radialMTFRec=radialMTFarrayRec(:,1)./radialMTFarrayRec(:,2);
radialMTFRec=cat(1,1,radialMTFRec);
axisradMTFRec=0:sqrt(2)*max(abs(axisMTF))/numbins:sqrt(2)*max(abs(axisMTF));

for i=1:sizMTFRes
   for j=1:sizMTFRes
      xposRes=i-peakxRes;
      yposRes=j-peakyRes;
      radiusRes=ceil(numbins*sqrt(xposRes^2+yposRes^2)/maxradiusRes);
      if (radiusRes==0)
      else
         radialMTFarrayRes(radiusRes,1)=radialMTFarrayRes(radiusRes,1)+MTFRes(i,j);
         radialMTFarrayRes(radiusRes,2)=radialMTFarrayRes(radiusRes,2)+1;
      end
   end
end


radialMTFRes=radialMTFarrayRes(:,1)./radialMTFarrayRes(:,2);
radialMTFRes=cat(1,1,radialMTFRes);
axisradMTFRes=0:sqrt(2)*max(abs(axisMTF))/numbins:sqrt(2)*max(abs(axisMTF));



figure;
set(gcf,'color','w');
plot(axisradMTFOrig,radialMTFOrig,axisradMTFRec,radialMTFRec,axisradMTFRes,radialMTFRes);
title('radial average MTF');
xlabel('c/deg');
ylabel('Radial Modulation Transfert Function');
legend('Original','recuperada','Residuo')

drawnow();
