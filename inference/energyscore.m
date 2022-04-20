function [score] = energyscore(y,mu,sig,k)
% ENERGYSCORE computes the multivariate energy score (as defined in
% Gneiting et al. (2008)) for predictions X on data Y, where X are K
% samples from the multivariate normal distribution defined by MU and
% SIG, e.g., the normal distribution defined by the ARs mean and std.
%
% The shape of Y, mu, and sig is [m,p], where m is the number of data
% points in time, and p the number of variables.
%

    % R. Marin 2022-04-20


    % reshape to proper size
    [m,p] = size(y);

    y = reshape(y,[m,1,p]);
    mu = reshape(mu,[m,1,p]);
    sig = reshape(sig,[m,1,p]);
    % sample the normal distribution
    X = mu + sig.*randn(m,k,p);

    % compute the Energy score
    norm1 = @(a,b) sqrt(nansum( (a-b).^2,1));
    norm2 = @(a) sqrt(nansum( diff(a,1,2).^2,1));
    score = squeeze(1/k*sum(norm1(X,y)) - 1/(2*(k-1))*sum(norm2(X)))';
end
