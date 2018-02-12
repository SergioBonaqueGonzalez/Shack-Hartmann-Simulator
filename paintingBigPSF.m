function [PSF,idealPSF]=paintingBigPSF(PSF,ML,SH,nFot,count,flag,PhotonNoise,ReadNoise,Exposuretime)
%{  
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This function paints the PSF of each microlense.
%} 
difference=length(PSF{1})-2*ML.radiusPixels;
BackGround=zeros(length(ML.AmplitudeMask)+difference);

%Joinning all PSF in the same figure
for i=1:length(PSF)
        BackGround(ML.coor(i,1):ML.coor(i,2)+difference, ML.coor(i,3):ML.coor(i,4)+difference) =...
            BackGround(ML.coor(i,1):ML.coor(i,2)+difference, ML.coor(i,3):ML.coor(i,4)+difference) + PSF{i};
end

if SH.paint==1 && count ==0
    figure
    subplot(1,2,1)
    imshow(BackGround,[])
    title('Image in the CCD (doubles spots are included)')
    set(gcf,'color','w');
    xlabel('pixels')
    ylabel('pixels')
end

BackGround(1:ceil(difference/2),:)=[];
BackGround(end-ceil(difference/2)+1:end,:)=[];
BackGround(:,1:ceil(difference/2))=[];
BackGround(:,end-ceil(difference/2)+1:end)=[];

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

if ML.SharedArea==1
    %Each area behind each microlens is takin into account for centroid
    %calculations. No matter if there exist energy from other microlens
    for i=1:length(PSF)
        idealPSF=PSF;
        PSF{i}=BackGround(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4));
    end
else
    idealPSF=[];
end

if SH.paint==1 && count==0
    subplot(1,2,2)
    imshow(BackGround.*ML.AmplitudeMask4Paint,[])
    title('Image in the CCD (the area of each microlenses is respected). The red point is the reference centroid')
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
