function [exp_var,exp_sample] = get_var(p,p0,N)
    %
    %returns the expected variance and number of evaluations
    %if p is used for allocation and p0 are the true probabilities
    %
    S = round(sum(1/p .* N));
    r = round(alloc(p,N));
    r(r==0)=1; 
    exp_var = sum(dilog(p0) ./ r) / S;
    exp_sd = sqrt(exp_var);

    exp_sample = sum(r ./ p0);
    

end 