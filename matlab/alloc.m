function r= alloc(p,N)
    r_n = sqrt(dilog(p) .* p);
    r_d = sqrt(dilog(p) ./ p);
    S = round(sum(1./p * N));
    r = S * r_n / sum(r_d);
end 