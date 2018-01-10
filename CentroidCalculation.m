function [cx, cy] = CentroidCalculation(image,method, display)
%{
Created by Sergio Bonaque-Gonzalez. Optical Engineer.
sergio.bonaque@um.es
This scripts calculates the centroid of an image by several ways.

INPUTS:
    Imagen= problem image. Any size is valid. 
    Manera= flag that indicates how the centroid will be calculated.
    display = 1, flag which indicates wheter paint or not the image and the centroid, pausing for a while.

OUTPUTS:
    [cx,cy]= coordinates of the centroid.
%}
cx=0;
cy=0;

%Only positive matrix are accepted as an input image
if min(min(image)) < 0
    disp('This function does not accept negative values in the input image');
    return
end

[size_y,size_x]=size(image);

%The special case where all values are zero is resolved here: 
if sum(sum(image))==0 
     cx=length(sum(image))/2;
     cy=length(sum(image,2)')/2;


% Method 1: classical way.
elseif method == 1    
    cx=(sum(image)*(1:length(sum(image)))')/sum(sum(image));
    cy=(sum(image,2)'*(1:length(sum(image,2)'))')/sum(sum(image));


%Method 2.
elseif method ==2
    % First,the maximum is calculated
    [~, indd] = max2(image);
    maximo_x=indd(2); 
    maximo_y=indd(1);
    vector_x= -maximo_x+1:(size_x-maximo_x);
    vector_y= -maximo_y+1:(size_y-maximo_y);
    delta_x=(sum(image)*(vector_x'))/sum(sum(image));
    delta_y=(sum(image,2)'*(vector_y)')/sum(sum(image));    
    %Finally, displacements from the maximum are calculated.
    cx= maximo_x+delta_x;
    cy= maximo_y+delta_y;
    
    
%Method 3: normalization of the peak
elseif method ==3
    [vall, ~] = max2(image);  
    imagen2=image./vall;
    cx=(sum(imagen2)*(1:length(sum(imagen2)))')/sum(sum(imagen2));
    cy=(sum(imagen2,2)'*(1:length(sum(imagen2,2)'))')/sum(sum(imagen2));  
    
    
%Method 4: setting a threshold in the average of the border, plus 3 times the standard deviation. After that, the centroid is calculated.  
elseif method ==4
    marco=[image(1,:) image(:,size_x)' image(size_y,:) image(:,1)'];
    media_marco=mean(marco);
    std_marco=std(marco);
    umbral = media_marco+3*std_marco;
    imagen2=image-umbral;
    imagen2(imagen2<0)=0; 
    cx=(sum(imagen2)*(1:length(sum(imagen2)))')/sum(sum(imagen2));
    cy=(sum(imagen2')*(1:length(sum(imagen2')))')/sum(sum(imagen2'));   
    
    
%Method 5: setting a threshold in the average of the whole image. Useful when borders have no information 
elseif method ==5
    roi=image(image > 0); 
    umbral = mean(roi(:));
    imagen2=image-umbral;
    imagen2(imagen2<0)=0;  
    cx=(sum(imagen2)*(1:length(sum(imagen2)))')/sum(sum(imagen2));
    cy=(sum(imagen2,2)*(1:length(sum(imagen2,2)'))')/sum(sum(imagen2));        


else
 disp('Not defined method for calculate the centroid');  
 return
end


%Painting
if method ==4
    if display == 1
        subplot(1,2,1);imshow(image, []);
        line([cx-1,cx+1],[cy ,cy],'color','b'); 
        line([cx,cx],[cy-1 ,cy+1],'color','b'); 

        if method ==4
            subplot(1,2,2);
            plot(image);
            hold;plot(umbral*ones(size(image)),'o'); 
        end
        shg;
        zoom(4);
        pause(0.3);
        close;
    end  
else
    if display == 1
        imshow(image, []);
        line([cx-1,cx+1],[cy ,cy],'color','b'); 
        line([cx,cx],[cy-1 ,cy+1],'color','b'); 
        shg;
        zoom(4);
        pause(0.1);
    end

end  
end

    

