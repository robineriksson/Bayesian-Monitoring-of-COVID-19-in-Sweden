function beta = getC19beta(rates,R,interp)
%GETC19BETA Compute beta for Covid-19 model.
%   beta = GETC19BETA(rates,R,interp) computes the value of beta given rate
%   parameters rates(n) and possibly time varying R. Daily time
%   resolution of the rates is required if R varies daily. interp
%   offers two different intepretations of the R0 defnintion, 2 being
%   suggested
%
%.  interp = 1 give the sqrt(R0) intepretation: the environmental
%   compartment is considered an infecitous state, or phi is the
%   adjusted infected levels. While interp = 2 gives the R0
%   intepretation: E is free infectious levels, bacteria in the
%   environment.
%
%   See also GETC19R0.

% H. Runvik 2021-02-22

if interp ~= 1 & interp ~= 2
    error('currently we only support interpretation 1 or 2');
end

beta_denom = rates.thetaE_./rates.sigma+ ...
           (rates.F0ave+(1-rates.F0ave).*rates.F1ave)./rates.gammaI+ ...
           (1-rates.F0ave).*rates.thetaA_./rates.gammaA;

if interp == 1
    beta = R.^2 ./ beta_denom;
elseif interp == 2
    beta = R ./ beta_denom;
end
