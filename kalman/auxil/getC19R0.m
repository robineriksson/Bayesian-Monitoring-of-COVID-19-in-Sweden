function R0 = getC19R0(rates,D, interp)
%GETC19R0 Compute R0 for Covid-19 model.
%   R0 = GETC19R0(rates,D,interp) computes the spatial effective value
%   of R0 given rate parameters rates(n) and transmission operator
%   D. interp offers two different intepretations of the R0
%   defnintion, 2 being suggested.
%
%   interp = 1 give the sqrt(R0) intepretation: the environmental
%   compartment is considered an infecitous state, or phi is the
%   adjusted infected levels. While interp = 2 gives the R0
%   intepretation: E is free infectious levels, bacteria in the
%   environment.
%
%   D can be left empty '[]', if only considering the direct conversion.
%
%   Example:
%     rates = priorenger(10,[3 3]);
%     load Dcounties
%     R0 = getC19R0(rates,D,1);
%
%     % compare:
%     [R0; rates.R0]
%
%   See also GETC19SYST.

% S. Engblom 2023-03-03 (Turns out there are two interpretations!)
% S. Engblom 2020-10-12 (Settled on formula)
% S. Engblom 2020-10-07

% note: does not depend on D nor on the time-discrete formulation:

if interp ~= 1 & interp ~= 2
    error('currently we only support interpretation 1 or 2');
end



R0 = rates.beta.* ...
     (rates.thetaE_./rates.sigma+ ...
     (rates.F0ave+(1-rates.F0ave).*rates.F1ave)./rates.gammaI+ ...
     (1-rates.F0ave).*rates.thetaA_./rates.gammaA);

if interp == 1
    R0 = sqrt(R0);
end

return

% see below for the derivation.

% $$$ % rates
% $$$ beta = rates.beta(n);
% $$$ sigma = rates.sigma(n);
% $$$ gammaI = rates.gammaI(n);
% $$$ gammaA = rates.gammaA(n);
% $$$ thetaI = rates.thetaI(n);
% $$$ thetaA = rates.thetaA(n);
% $$$ thetaE = rates.thetaE(n);
% $$$ rho = rates.rho(n);
% $$$ F0 = rates.F0ave(n);
% $$$ F1 = rates.F1ave(n);
% $$$
% $$$ % (thetaI_ = 1 implicitly understood)
% $$$ thetaA_ = thetaA/thetaI;
% $$$ thetaE_ = thetaE/thetaI;
% $$$
% $$$ % this scaling uses that thetaI = rho:
% $$$ erho = -expm1(-rho); % (skip division by rho here)
% $$$
% $$$ % form the Next Generation Operator NG: for states x = [E A I phi],
% $$$ % i.e., x' = (F-V)x, or x(k+1) = (I+F-V)x(k), in terms of which NG =
% $$$ % F\V
% $$$ F = [0 0 0 beta; ...
% $$$      zeros(2,4); ...
% $$$      erho*thetaE_ erho*thetaA_ erho 0];
% $$$ V = [sigma 0 0 0; ...
% $$$      -sigma*(1-F0) gammaA 0 0; ...
% $$$      -sigma*F0 -gammaA*F1 gammaI 0; ...
% $$$      0 0 0 erho];
% $$$ % (see https://en.wikipedia.org/wiki/Next-generation_matrix)
% $$$ R0 = max(abs(eig(F/V)));

% symbolic check:
syms beta erho thetaE_ thetaA_
syms sigma gammaA gammaI
syms F0 F1

% x = [E A I phi]
F = [0 0 0 beta; ...
     zeros(2,4); ...
     erho*thetaE_ erho*thetaA_ erho 0];
V =  [sigma 0 0 0; ...
     -sigma*(1-F0) gammaA 0 0; ...
     -sigma*F0 -gammaA*F1 gammaI 0; ...
     0 0 0 erho];
lam = eig(F/V);

% result after some manual polishing:
sqrt(beta* ...
     (thetaE_/sigma + ...
      (F0+(1-F0)*F1)/gammaI + ...
      thetaA_*(1-F0)/gammaA));
res = lam(4)^2-beta* ...
      (thetaE_/sigma + ...
       (F0+(1-F0)*F1)/gammaI + ...
       thetaA_*(1-F0)/gammaA);
simplify(res) % = 0

% move source terms of phi to V - so E alone is the state-at-infection
F = [0 0 0 beta; ...
     zeros(3,4)];
V =  [sigma 0 0 0; ...
     -sigma*(1-F0) gammaA 0 0; ...
     -sigma*F0 -gammaA*F1 gammaI 0; ...
      [0 0 0 erho]- ...
      [erho*thetaE_ erho*thetaA_ erho 0]];
lam = eig(F/V); % = lam(4)^2, so square of the above
res = lam(4)-beta* ...
      (thetaE_/sigma + ...
       (F0+(1-F0)*F1)/gammaI + ...
       thetaA_*(1-F0)/gammaA);
simplify(res) % = 0

% move source term of E to V - so phi alone is the effective number of
% infected individuals - the state-at-infection
F = [0 0 0 0; ...
     zeros(2,4); ...
     erho*thetaE_ erho*thetaA_ erho 0];
V =  [sigma 0 0 -beta; ...
     -sigma*(1-F0) gammaA 0 0; ...
     -sigma*F0 -gammaA*F1 gammaI 0; ...
     0 0 0 erho];
lam = eig(F/V); % = lam(4)^2, so square of the above
res = lam(4)-beta* ...
      (thetaE_/sigma + ...
       (F0+(1-F0)*F1)/gammaI + ...
       thetaA_*(1-F0)/gammaA);
simplify(res) % = 0

% (not needed as a conservative D does not affect R0:)

% $$$ % effective transmission matrix
% $$$ Deff = D-diag(sum(D,1));
% $$$
% $$$ % "networkify" the model
% $$$ nstate = 4;
% $$$ T = zeros(nstate); T(4,4) = 1; % (phi)
% $$$ D = kron(rates.lambda(n)*Deff,T);
% $$$
% $$$ Id_lan = speye(size(Deff));
% $$$ F = kron(Id_lan,F);
% $$$ V = kron(Id_lan,V)-D;
% $$$ R0 = max(abs(eigs(F/V)));
