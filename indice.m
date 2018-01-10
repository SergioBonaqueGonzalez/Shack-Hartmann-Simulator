function [ngra mazi] = indice(jfin)
%Created by Sergio Bonaque-Gonzalez. Optical Engineer.
%   sergio.bonaque@um.es
% This is an auxiliary function for compute Zernike coefficients.

    n_gra = zeros(1, jfin + 10);
    m_azi = zeros(1, jfin + 10);

    % RADIAL DEGREES AND AZIMUTHAL FRENCUENCY ASSOCIATED TO POLYNOMIAL J INDEX.

    n = 1;
    j = 0.0;
    while (j <= jfin)
        for mm = 0:n
            j = ((n * (n + 1.0)) / 2.0) + mm + 1.0;
            if (paridad(n) ~= paridad(mm)) 
                m = mm + 1;
            else
                m = mm;
            end
            
            n_gra(j) = n;
            m_azi(j) = m;
        end
        
        n = n + 1;
    end
        
    ngra = n_gra(jfin);
    mazi = m_azi(jfin);
end
