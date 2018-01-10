function p = paridad(x)
% Auxiliary function for calculating Zernike coefficients. 
    if (mod(x, 2) == 0) 
        p = true;
    else 
        p = false; 
    end
end
