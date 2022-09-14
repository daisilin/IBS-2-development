function v_ratio= compute_vratio_s16(p)
    v_ratio=    (sum(sqrt(dilog(p)./p)))^2 * (1./sum(dilog(p))) * (1./ sum(1./p));
end 
