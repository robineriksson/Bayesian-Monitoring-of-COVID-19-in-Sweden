function [xcov,xmean]=covRec(x,oldcov,oldmean,w)
%COVREC update the covariance matrix by recursion, when a new samples
%   arrives.
%   [XCOV,XMEAN]=COVREC(X,OLDCOV,OLDMEAN,W) returns the recursive mean and
%   the recursive sample covariance, using the new sample X and the old covariance
%   OLDCOV and the old mean OLDMEAN, on the W:th sample.
%
%   X       - a [1 x n] matrix. 
%   OLDCOV  - previous step sample covariance.
%   OLDMEAN - previous step sample mean.
%   W       - sample number, e.g., w \in {1:n}
%
%   initial step, i.e., oldcov = oldmean = [] will return a zero
%   initalization.


% R. Eriksson 2021-01-25


[p, n] = size(x);
if p > 1
  error('x must be a column matrix');
end

if w > 0 && ~isempty(oldcov) % update recursively
  xmean  = oldmean + 1/(w)*(x-oldmean);
  
  % biased estimate
  %xcov = (w-1)/w*oldcov + (w-1)/w^2*(x - oldmean)'*(x-oldmean);
  
  % unbiased estimate
  xcov =  oldcov + 1/(w-1)*( (w-1)/(w)*((x-oldmean)'*(x-oldmean)) - oldcov);
else % initialize
  xmean = x;
  xcov = zeros(n);
end


end