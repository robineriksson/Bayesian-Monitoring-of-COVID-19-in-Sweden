function R0 = getC19R0(rates,D)
%GETC19R0 Compute R0 for Covid-19 model.
%   R0 = GETC19R0(rates,D) computes the spatial effective value of R0
%   given rate parameters rates(n) and transmission operator D.
%
%   Example:
%     rates = priorenger(10,[3 3]);
%     load Dcounties
%      R0 = getC19R0(rates,D);
%
%     % compare:
%     [R0; rates.R0]

%   See also GETC19SYST.

% S. Engblom 2020-10-12 (Settled on formula)
% S. Engblom 2020-10-07

% note: does not depend on D nor on the time-discrete formulation:
R0 = sqrt(rates.beta.* ...
          (rates.thetaE_./rates.sigma+ ...
           (rates.F0ave+(1-rates.F0ave).*rates.F1ave)./rates.gammaI+ ...
           (1-rates.F0ave).*rates.thetaA_./rates.gammaA));

return;


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
% $$$ 
% $$$ % symbolic check:
% $$$ beta = sym('beta');
% $$$ erho = sym('erho');
% $$$ thetaE_ = sym('thetaE_');
% $$$ thetaA_ = sym('thetaA_');
% $$$ F = [0 0 0 beta; ...
% $$$      zeros(2,4); ...
% $$$      erho*thetaE_ erho*thetaA_ erho 0];
% $$$ sigma = sym('sigma');
% $$$ gammaA = sym('gammaA');
% $$$ gammaI = sym('gammaI');
% $$$ F0 = sym('F0');
% $$$ F1 = sym('F1');
% $$$ V =  [sigma 0 0 0; ...
% $$$      -sigma*(1-F0) gammaA 0 0; ...
% $$$      -sigma*F0 -gammaA*F1 gammaI 0; ...
% $$$      0 0 0 erho];
% $$$ lam = eig(F/V);
% $$$ % result after some manual polishing:
% $$$ sqrt(beta* ...
% $$$      (thetaE_/sigma + ...
% $$$       (F0+(1-F0)*F1)/gammaI + ...
% $$$       thetaA_*(1-F0)/gammaA));
% $$$ res = lam(4)^2-beta* ...
% $$$       (thetaE_/sigma + ...
% $$$        (F0+(1-F0)*F1)/gammaI + ...
% $$$        thetaA_*(1-F0)/gammaA);
% $$$ simplify(res) % = 0

% *** not needed as a conservative D does not affect R0:

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
