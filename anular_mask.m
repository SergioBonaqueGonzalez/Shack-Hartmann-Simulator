function mask = anular_mask(dimx, dimy, X_cent, Y_cent, R_int, R_ext)
%Created by Juan Valdivia and modified by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This function calculates a anular and circular masks.

    mask  = zeros(dimx, dimy);
    dimx_ = dimx - 1;
    dimy_ = dimy - 1;
    
    for p = 0:dimy_
        for i = 0:dimx_
            radio = realsqrt(((i - X_cent) * (i - X_cent)) + ((p - Y_cent) * (p - Y_cent)));

            if ((radio >= R_int) && (radio <= R_ext))
                mask(p + 1, i + 1) = 1.0;
            end
        end
    end
end

	
