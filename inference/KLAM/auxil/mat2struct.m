function rates = mat2struct(ratesMat,ratenames,hyp,fix,nslab)
%RATES = mat2struct(RATESMAT,RATENAMES,LOHI,FIX,TOEXP,NSLAB,HYP) From
%   a matrix with some selected rates, generate a struct, and the missing
%   parameters
%
%   RATESMAT  = the matrix which hold the parameters we want to store in the
%               struct
%   RATESNAME = the name of the rates in the matrix
%   FIX       = true if, the 'fixed' parameter space is used or not
%   NSLAB     = the number of slabs used
%   HYP       = hyperparameters loaded from preceeding prior.

% R. Eriksson 210522 (modernized into to couple with priorenger)
% R. Eriksson 201006

switch nargin
  case 3
    fix = true;
    nslab = 1;
  case 4
    nslab = 1;
  case 5
    % nothing
  otherwise, error('Incorrect number of arguments.');
end


Nreplicas = size(ratesMat,2);
rates = struct();

% we assume that there are slabdependencies
uniqueNames = fsetop('unique',ratenames);
for i = 1:numel(uniqueNames)
 name = uniqueNames{i};
 rate = ratesMat(strcmp(ratenames,name),1:Nreplicas);
 rates.(name) = rate;
end


if fix
 rates.gammaI = rates.sigma;
end
rates.gammaA = rates.gammaI; % always fixed.

% add bacterial decay
rates.rho = -log(0.5)./rates.half_life;
rates.A2I = zeros(size(rates.E2I));

try % check if lambda is included or not in the matrix.
 rates.lambda;
catch
 rates.lambda = hyp.lambdaPQ(1)+diff(hyp.lambdaPQ)*...
   betarnd(hyp.lambdaHyp(1),hyp.lambdaHyp(2),nslab,Nreplicas);
end

try % check if lambda is included or not in the matrix.
 rates.k_sew;
catch
    rates.k_sew = NaN(nslab,Nreplicas); %*** maybe switch to NaN
end

% add shedding
rates.thetaI_ = ones(1, Nreplicas);
rates.thetaI = rates.thetaI_.*rates.rho;
rates.thetaA = rates.thetaA_.*rates.thetaI;
rates.thetaE = rates.thetaE_.*rates.thetaI;

% SIR_MORT, HOSP_MORT, F3scale: known and fixed values
SIR_MORT = hyp.SIR_MORT;
% HOSP_MORT
HOSP_MORT = hyp.HOSP_MORT;

% Demographic averages:
% F0ave, F1ave, F2ave, F2dave, F3ave, F3dave, F4ave
rates.F0ave = repmat(rates.E2I,nslab,1);

rates.F1ave = repmat(rates.A2I,nslab,1);

rates.F4ave = repmat(SIR_MORT,nslab,Nreplicas);

rates.F3ave = repmat(rates.IC_HOSP,nslab,1);

rates.F3dave = repmat(SIR_MORT*HOSP_MORT,nslab,Nreplicas);

rates.F2dave = max(rates.IFR./rates.E2I-...
  (HOSP_MORT+rates.IC_HOSP)*SIR_MORT.*rates.HOSP./...
  (1-rates.IC_HOSP*(1-SIR_MORT)),0);

rates.F2ave = repmat(rates.HOSP,nslab,1);


% deduced contraction rate beta
rates.beta = getC19beta(rates,rates.R0,hyp.interp);

end
