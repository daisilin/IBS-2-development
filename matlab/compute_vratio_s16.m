function v_ratio= compute_vratio_s16(p)
    v_ratio= (sum(sqrt(dilog(p)./p)))^2 * (1./sum(dilog(p))) * (1./ sum(1./p));
%     f1 = 0;
%     f2 = 0;
%     f3 = 0;
%     for i = 1:length(p)
%         p_i = p(i);
%         f1 = f1+ sqrt(dilog(p_i)/p_i);
%         f2 = f2+ dilog(p_i);
%         f3 = f3 + 1/p_i;
%     end 
%     v_ratio = f1^2 * 1/f2 * 1/f3;
end 
