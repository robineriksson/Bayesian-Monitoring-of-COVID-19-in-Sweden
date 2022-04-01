function beta = getC19beta(rates,R)
%GETC19BETA Compute beta for Covid-19 model.
%   beta = GETC19BETA(rates,R) computes the value of beta given rate
%   parameters rates(n) and possibly time varying R. Daily time resolution
%   of the rates is required if R varies daily.
%   
%   See also GETC19R0.

% H. Runvik 2021-02-22

beta = R.^2./(rates.thetaE_./rates.sigma+ ...
           (rates.F0ave+(1-rates.F0ave).*rates.F1ave)./rates.gammaI+ ...
           (1-rates.F0ave).*rates.thetaA_./rates.gammaA);


