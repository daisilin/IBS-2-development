function [L,L_trials]=compute_nll_psycho(stim,resp,theta)
    sigma=exp(theta(1));
    bias=theta(2);
    lambda=theta(3);
    L_trials=log(lambda/2+(1-lambda)*((resp==1).*normcdf(-(stim-bias)/sigma)+(resp==0).*normcdf((stim-bias)/sigma)));
    L = -sum(L_trials);
end