function Nreps_dynamic=Nreps_dynamic_estimate(p_current, Nsamples)
S_budget = round(sum(1./p_current * Nsamples));%estimate the total number of sample (budget)
Nreps_dynamic  = round(S_budget * (1./sum(sqrt(dilog(p_current)./p_current)))*sqrt(p_current.*dilog(p_current)));

end 