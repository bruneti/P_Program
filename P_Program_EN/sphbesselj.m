function jl = sphbesselj(l, x)
    jl = sqrt(pi./(2*x)) .* besselj(l + 0.5, x);
    % handle x = 0 to avoid division by zero
    jl(x==0) = double(l==0);
end
