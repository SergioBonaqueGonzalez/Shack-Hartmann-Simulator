function zz = zernike(j, dim)

    temp = zeros(dim, dim);
    xx   = zeros(dim, dim);
    yy   = zeros(dim, dim);
    tdim = dim / 2.;
    
    for i = 1:dim
        for p = 1:dim
            xx(i, p) = p - 1 - tdim;
            yy(i, p) = i - 1 - tdim;
        end
    end

    tdim   = (dim / 2) - 1;
    ro     = realsqrt((xx .* xx) + (yy .* yy)) ./ tdim;
    teta   = atan2(yy, xx);
    [n, m] = indice(j);

%     mask = anular_mask(dim, dim, dim / 2, dim / 2, 0, (dim / 2) - 1);
    mask = mascaraCircular(1,dim);

    if (m == 0)
        temp = realsqrt(n + 1) .* zer_rad(ro, n, m) .* mask;
    end

    if (m ~= 0)
        if (paridad(j) == true)
            temp = realsqrt((2 * n) + 2) .* zer_rad(ro, n, m) .* cos(m .* teta) .* mask;
        end

        if (paridad(j) == false)
            temp = realsqrt((2 * n) + 2) .* zer_rad(ro, n, m) .* sin(m .* teta) .* mask;
        end
    end

    zz = temp;