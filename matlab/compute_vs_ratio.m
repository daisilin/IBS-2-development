function vs_ratio= compute_vs_ratio(p0,p_vec)
        
vs_ratio = (sum(sqrt(dilog(p0).*p0).*dilog(p0)) * sum(sqrt(dilog(p_vec).*p_vec)./p0)) / (sum(sqrt(dilog(p_vec).*p_vec).*dilog(p0)) * sum(sqrt(dilog(p0).*p0)./p0)); 
end 
