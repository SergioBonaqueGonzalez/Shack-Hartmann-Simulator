function outputImage = lensletSimulation(S,L,lambda,focalLength,pixelSize)

% Created by Juan Trujillo.
% jtrujils@ull.es
%	This simulates individual lenslet by Fraunhofer propagation
% 
% INPUTS:
%     S = input complex field;
%     L = side length of input complex field
%     lambda = wave length
%     focalLength = Lens focal length in meters 
%     pixelSize = sensor pixel size in meters
% 
% OUTPUTS:
%     outputImage = Output image at sensor plane


    [propagated, L2] =FraunhoferPropagation(S,L,lambda,focalLength); %L2 represent the real size of the PSF in meters
    outputImage = (abs(propagated).^2);
    outputImage = imresize(outputImage,(L2/length(propagated))/pixelSize,'bilinear'); %Scaling the pixel size of the propagated field to the correct one.
    remainingPixels = (size(outputImage,1) -  size(S,1))/2;

    if remainingPixels < 0 %padding
        if mod(remainingPixels,2) ~= 0
            outputImage = imtranslate(outputImage,[0.5,0.5]);
            remainingPixels = remainingPixels - 0.5;                
        end
        outputImage = padarray(outputImage,abs([remainingPixels remainingPixels]),0,'both');
%     elseif remainingPixels > 0 %crop
%         if mod(remainingPixels,2) ~= 0
%             outputImage = imtranslate(outputImage,[0.5,0.5]);
%             remainingPixels = remainingPixels - 0.5;
%         end
%         outputImage = imcrop(outputImage,[remainingPixels+1,remainingPixels+1,size(S,1)-1,size(S,1)-1]);
    end

end