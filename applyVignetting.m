function outputImage = applyVignetting(focalLength,pixelSize,microlensPixels,microlensThickness,image)


    pupilDiameter = pixelSize*microlensPixels;

    alphaMax = atan(0.5*pupilDiameter/focalLength);

    [u,v] = meshgrid(linspace(-alphaMax,alphaMax,size(image,1)),linspace(-alphaMax,alphaMax,size(image,2)));
    u = abs(u);
    v = abs(v);


    effectivePupilArea = zeros(size(image));
    for i = 1:size(image,1)
        for j = 1:size(image,1)
           currentAngle = sqrt(u(i,j)^2+v(i,j)^2);
           effectivePupilArea(i,j) =  1 - microlensThickness*tan(currentAngle);
        end
    end

    outputImage = image .* effectivePupilArea;