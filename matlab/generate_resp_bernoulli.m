function resp=generate_resp_bernoulli(stim,p,ind)
%GENERATE_RESP_BERNOULLI Generate responses for bernoulli.

if nargin < 3 || isempty(ind); ind = (1:size(stim,1))'; end

stim=stim(ind);
resp=rand(size(stim))<(p(ind));
    
end