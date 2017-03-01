function [cx, cy] = CalCentroides( imagen, manera, display)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
%This scripts calculates the centroid of an image by several ways.

%INPUTS:
%Imagen= problem image. Any size is valid. 
%Manera= flag that indicates how the centroid will be calculated.
%display = 1, flag which indicates wheter paint or not the image and the centroid, pausing for a while.

%OUTPUTS:
%[cx,cy]= coordinates of the centroid.

%SCRIPT
cx=0;
cy=0;

%Only positive matrix are accepted as an input image
if min(min(imagen)) < 0
    disp('This function does not accept negative values in the input image');
    return
end

[size_y,size_x]=size(imagen);

%The special case where all values are zero is resolved here: 
if sum(sum(imagen))==0 
     cx=length(sum(imagen))/2;
     cy=length(sum(imagen'))/2;


% Method 1: classical way.
elseif manera == 1    
    cx=(sum(imagen)*(1:length(sum(imagen)))')/sum(sum(imagen));
    cy=(sum(imagen')*(1:length(sum(imagen')))')/sum(sum(imagen'));


%Method 2.
elseif manera ==2
    % First,the maximum is calculated
    [vall, indd] = max2(imagen);
    maximo_x=indd(2); 
    maximo_y=indd(1);
    vector_x= -maximo_x+1:(size_x-maximo_x);
    vector_y= -maximo_y+1:(size_y-maximo_y);
    delta_x=(sum(imagen)*(vector_x'))/sum(sum(imagen));
    delta_y=(sum(imagen')*(vector_y)')/sum(sum(imagen'));    
    %Finally, displacements from the maximum are calculated.
    cx= maximo_x+delta_x;
    cy= maximo_y+delta_y;
    
    
%Method 3: normalization of the peak
elseif manera ==3
    [vall, indd] = max2(imagen);  
    imagen2=imagen./vall;
    cx=(sum(imagen2)*(1:length(sum(imagen2)))')/sum(sum(imagen2));
    cy=(sum(imagen2')*(1:length(sum(imagen2')))')/sum(sum(imagen2'));  
    
    
%Method 4: setting a threshold in the average of the border, plus 3 times the standard deviation. After that, the centroid is calculated.  
elseif manera ==4
    marco=[imagen(1,:) imagen(:,size_x)' imagen(size_y,:) imagen(:,1)'];
    media_marco=mean(marco);
    std_marco=std(marco);
    umbral = media_marco+3*std_marco;
    imagen2=imagen-umbral;
    imagen2(imagen2<0)=0; 
    cx=(sum(imagen2)*(1:length(sum(imagen2)))')/sum(sum(imagen2));
    cy=(sum(imagen2')*(1:length(sum(imagen2')))')/sum(sum(imagen2'));   
    
    
%Method 5: setting a threshold in the average of the whole image. Useful when borders have no information 
elseif manera ==5
    roi=imagen(imagen > 0); 
    umbral = mean(roi(:));
    imagen2=imagen-umbral;
    imagen2(imagen2<0)=0;  
    cx=(sum(imagen2)*(1:length(sum(imagen2)))')/sum(sum(imagen2));
    cy=(sum(imagen2')*(1:length(sum(imagen2')))')/sum(sum(imagen2'));        


else
 disp('Not defined method for calculate the centroid');  
 return
end


%Painting
if manera ==4
    if display == 1
        figure;
        subplot(1,2,1);imshow(imagen, []);
        line([cx-1,cx+1],[cy ,cy],'color','b'); 
        line([cx,cx],[cy-1 ,cy+1],'color','b'); 

        if manera ==4
            subplot(1,2,2);
            plot(imagen);
            hold;plot(umbral*ones(size(imagen)),'o'); 
        end
        shg;
        zoom(4);
        pause(0.3);
        close;
    end  
else
    if display == 1
        imshow(imagen, []);
        line([cx-1,cx+1],[cy ,cy],'color','b'); 
        line([cx,cx],[cy-1 ,cy+1],'color','b'); 
        shg;
        zoom(4);
        pause(0.1);
    end

end  
end

    

