function Nreps = estimate_Nreps(Nsamples,p_current,lookup_logp,dilog_p)
    S_budget = round(sum(1./p_current * Nsamples));%estimate the total number of sample (budget)
    dilog_lookup = interp1(lookup_logp,dilog_p,log(p_current));
    Nreps = round(S_budget * (1./sum(sqrt(dilog_lookup./p_current)))*sqrt(p_current.* dilog_lookup ));
    Nreps(Nreps==0) = 1;
end 