function jl = sphbesselj(l, x)
    jl = sqrt(pi./(2*x)) .* besselj(l + 0.5, x);
    % caso x = 0 para evitar la división por cero
    jl(x==0) = double(l==0);
end