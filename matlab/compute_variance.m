function nll_var_vec_i= compute_variance(K_active,ind_active,Nreps)
Ktab = -(psi(1,1:max(K_active(:)))' - psi(1,1));
% nll_var_vec_i = nansum(Ktab(K_active)./Nreps(ind_active).^2);
nll_var_vec_i = nansum(Ktab(K_active));

end 