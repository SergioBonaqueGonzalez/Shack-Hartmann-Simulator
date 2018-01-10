function distortedImage = applySeidelDistortion(focalLength,pixelSize,fieldDistortion,microlensPixels,image)

% Created by Juan Trujillo.
% jtrujils@ull.es
%	This function applies Seidel field distortion to an image produced by a
%	lens of known focal length, in a sensor with known pixel size.
% 
% INPUTS:
%     focalLength = Lens focal length in meters 
%     pixelSize = sensor pixel size in meters
%     fieldDistortion = Coefficient W311 of Seidel aberrations in meters
%     image = input image to be distorted
%     microlensPixels = number of linear pixel sin a microlens
% 
% OUTPUTS:
%     distortedImage = Output distorted image

    imageSpaceNormalizedWidth = 0.5*size(image,1)/microlensPixels;

    conversionFactor = (pixelSize^2)/focalLength;
    [X,Y] = meshgrid(linspace(-0.5,0.5,size(image,1)),linspace(-0.5,0.5,size(image,2)));
    [u,v] = meshgrid(linspace(-imageSpaceNormalizedWidth,imageSpaceNormalizedWidth,size(image,1)),...
        linspace(-imageSpaceNormalizedWidth,imageSpaceNormalizedWidth,size(image,2)));
    distortionOPD = zeros(size(image));
    for ii = 1:size(image,1)
        for jj = 1:size(image,2)
            distortionOPD(ii,jj) = seidel_5(u(ii,jj),v(ii,jj),X(ii,jj),Y(ii,jj),0,0,0,0,0,fieldDistortion);
        end
    end
    [xDistortion, yDistortion] = imgradientxy(distortionOPD/conversionFactor);
    d(:,:,1) = xDistortion;
    d(:,:,2) = yDistortion;
    distortedImage = imwarp(image, -d); 
    distortedImage = sum(image(:))*distortedImage/sum(distortedImage(:));

end


function[w]=seidel_5(u0,v0,X,Y,...
wd,w040,w131,w222,w220,w311)
    % seidel_5
    % Compute wavefront OPD for first 5 Seidel wavefront
    % aberration coefficients + defocus
    %
    %
    % u0,v0 - normalized image plane coordinate
    % X,Y - normalized pupil coordinate arrays
    % (like from meshgrid)
    % wd-defocus; w040-spherical; w131-coma;
    % w222-astigmatism; w220-field curvature;
    % w311-distortion

    beta=atan2(v0,u0); % image rotation angle
    u0r=sqrt(u0^2+v0^2); % image height 
    % Wavefront Aberrations 145
    % rotate grid
    Xr=X*cos(beta)+Y*sin(beta);
    Yr=-X*sin(beta)+Y*cos(beta);

    % Seidel polynomials
    rho2=Xr.^2+Yr.^2;
    w=wd*rho2+...
    w040*rho2.^2+...
    w131*u0r*rho2.*Xr+...
    w222*u0r^2*Xr.^2+...
    w220*u0r^2*rho2+...
    w311*u0r^3*Xr;

end