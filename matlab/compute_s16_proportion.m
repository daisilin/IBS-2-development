function v_ratio= compute_s16_proportion(p_vec,prop_vec)
    f1=0;
    f2=0;
    f3=0;
    for i = 1:length(p_vec)
        f1 = f1 + prop_vec(i)* (1/p_vec(i));
        f2 = f2 + prop_vec(i)* dilog(p_vec(i));
        f3 = f3 + prop_vec(i)* sqrt(dilog(p_vec(i))/p_vec(i));
    end 
    
    v_ratio= 1/f1 * 1/f2 * f3^2;
end 
