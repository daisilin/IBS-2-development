function [Nreps,S_budget] = estimate_Nreps(Nsamples,lookup_logp,dilog_p,p_vec,fix_budget,p_true)
    if  isscalar(p_vec)
        Nreps = Nsamples;
        S_budget = NaN;
    elseif fix_budget == true 
%         p_half = ones(length(p_vec),1)*0.5;
        S_budget = round(sum(1./p_true * Nsamples));%estimate the total number of sample (budget) using TRUE p
%         S_budget = 800;
        dilog_lookup = interp1(lookup_logp,dilog_p,log(p_true));
%         Nreps = round(S_budget * (1./sum(sqrt(dilog_lookup./p_true)))*sqrt(p_vec.* dilog_lookup ));
%         Nreps = round(S_budget * (1./sum(sqrt(dilog(p_true)./p_true))) * sqrt(p_vec.*dilog(p_vec))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs
        Nreps = round(S_budget * sqrt(p_vec.*dilog(p_vec)));
        scaling = S_budget/(sum(Nreps.*(1./p_true)));
        Nreps = round(Nreps*scaling);
        Nreps(Nreps==0) = 1;         
    else 
        S_budget = round(sum(1./p_vec * Nsamples));%estimate the total number of sample (budget)
        dilog_lookup = interp1(lookup_logp,dilog_p,log(p_vec));
        Nreps = round(S_budget * (1./sum(sqrt(dilog_lookup./p_vec)))*sqrt(p_vec.* dilog_lookup ));
        Nreps(Nreps==0) = 1;       
    end  
        
end 