function [PSF,idealPSF]=paintingPSF(PSF,ML,SH,nFot,count,flag,PhotonNoise,ReadNoise,Exposuretime)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function paints the PSF of each microlense.
As photon noise will not be used to calibration, flag indicates if this is
or not a calibration procedure 0=NO, 1=YES
%}
idealPSF=[];
BackGround=zeros(length(ML.AmplitudeMask));

%Joinning all PSF in the same figure
for i=1:length(PSF)
    BackGround(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4))=PSF{i};
end

%Introduction of Photon Noise
if PhotonNoise==1 && flag==1
    BackGround = poissrnd(nFot*BackGround/sum(BackGround(:)));
    for i=1:length(PSF)
        PSF{i}=BackGround(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4));
    end
end

%Introduction of Read Noise
if ReadNoise==1 && flag==1
    BackGround = readNoise(BackGround,SH,Exposuretime);
    for i=1:length(PSF)
        PSF{i}=BackGround(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4));
    end
end

if SH.paint==1 && count == 0
    figure
    subplot(1,2,1)
    imshow(BackGround,[])
    title('Image in the CCD')
    set(gcf,'color','w');
    xlabel('pixels')
    ylabel('pixels')
    
    subplot(1,2,2)
    imshow(BackGround.*ML.AmplitudeMask4Paint,[])
    title('Image in the CCD with centroid (red point) and asigned area for each microlens')
    xlabel('pixels')
    ylabel('pixels')
    hold on
    for i=1:length(PSF)
        rectangle('Position', [ML.coor(i,1) ML.coor(i,3) ML.radiusPixels*2 ML.radiusPixels*2],'LineWidth', 0.1, 'EdgeColor', 'b');
        %pause(0.01)
    end
    hold on
    plot(ML.coor(:,1)+ML.radiusPixels,ML.coor(:,3)+ML.radiusPixels,'o','MarkerSize',1,'MarkerEdgeColor','r')
    
    drawnow();
end
