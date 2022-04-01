function [rates,names,c19prior] = priorenger(X,fix,nslab,hypfile)
%PRIORENGER Prior parameter distribution/density for the EngEr Covid-19 model.
%
%   Case 1: prior sample (X is an integer or inf).
%
%   RATES = PRIORENGER(X,FIX,NSLAB,HYPFILE) returns a structure RATES with X
%   rate parameters according to FIX, for NSLAB number of slabs. The
%   parameters we assume change in slabs are beta (R0), IFR (infection
%   fatality rate), and lambda (transport), which means that these will be
%   returned as NSLAB-by-X matrices. HYP holds the hyperparameters for the
%   rate priors.
%
%   When X = Inf, the mean of the prior is rather returned.
%
%   Case 2: prior density (X is a struct).
%
%   RATES = PRIORENGER(X,[...]) returns a structure RATES with the
%   marginal density of the parameters in X, assuming X holds the same
%   structure as RATES. RATES.SUMLOG is the sum(log) of all marginal
%   densities. RATES.SUMLOG_DYNAMIC is the sum(log) per slab for the
%   dynamic parameters (see above). RATES.SUMLOG = sum(log({fixed
%   rates} \ {dynamic rates})) + RATES.SUMLOG_DYNAMIC.
%
%   FIX = true/{false}, a flag that if true reduces the parameter
%   space such that sigma = gammaI = gammaA, i.e., only one random
%   variable is sampled per group reducing the prior density
%   dimensionality by 2.
%
%   NSLAB = integer {1}, assumes dynamic rate parameters in NSLAB
%   slabs. See C19FILT for how the slabs are handeled.
%
%   HYPFILE = {'inference/c19prior'}, character leading pointing to file of
%   hyperparemeters see GENPRIOR.m, either left empty, and then a set of
%   standard values are used, otherwise, a file is loaded. For example:
%      % construct hyperparameter struct
%      genprior;
%      save('inference/c19prior','c19prior');
%
%   [RATES,NAMES] = PRIORENGER also returns rate names and values,
%   respectively.
%
%   [[...],C19PRIOR] = PRIORENGER, the third output variable, C19PRIOR are
%   the hyperparameters used. A case when this output is needed using
%   mat2struct at a later stage, see Example 3.
%
%   Example 1:
%     % Get a single prior sample, and its density
%     rates = priorenger(1);
%     dens = priorenger(rates);
%     sumloglike = dens.sumlog
%
%   Example 2:
%     % Get a 10 prior samples with 5 slabs, with a preloaded F
%     rates = priorenger(10,false,5);
%
%   Example 3:
%     % get a prior sample, keep the hyperparameters
%     file = 'inference/c19prior';
%     [rates,names,hyp] = priorenger(1);
%     % extract the vector that would be used in AM
%     [mat,ratenames,meta] = struct2mat(rates);
%     % convert back into a struct | remember the hyperparameters
%     rates_ = mat2struct(mat,ratenames,hyp,fix,nslab);
%     rates_.meta = meta;
%     isequal(rates,rates_) % is true.
%


% R. Eriksson 2020-05-07 (removed F, and scaleF call)
% R. Eriksson 2020-05-05 (removed lohi and chunks)
% R. Eriksson 2020-10-26 (introduced slab proposals)
% R. Eriksson (merge with priorengerDens)
% R. Eriksson 2020-09-22 (update of parameters distribution)
% S. Engblom 2020-04-17


switch nargin
  case 0
    X = inf;
    fix = false;
    nslab = 1;
    hypfile = mfilename('fullpath');
    hypfile = [hypfile(1:end-10) 'c19prior'];
  case 1
    fix = false;
    nslab = 1;
    hypfile = mfilename('fullpath');
    hypfile = [hypfile(1:end-10) 'c19prior'];
  case 2
    nslab = 1;
    hypfile = mfilename('fullpath');
    hypfile = [hypfile(1:end-10) 'c19prior'];
  case 3
    hypfile = mfilename('fullpath');
    hypfile = [hypfile(1:end-10) 'c19prior'];
  case 4
  otherwise, error('Incorrect number of arguments.');
end

% To reduce the number of lines, define Nreplicas to hold the number
% of replicas, and inf if we only need 1.
if isstruct(X) % probability
  Nreplicas = size(X.R0,2);
elseif X == inf % mean sample
  Nreplicas = 1;
else % sample
  Nreplicas = X;
end

% turn Matlab native pdf's on/off
nativepdf = false;
if nativepdf
  betadens = @(x,a,b) betapdf(x,a,b);
  unifdens = @(x,a,b) unifpdf(x,a,b);
  truncLogNormpdf = @(dist,x) pdf(dist,x);
else
  betadens = @(x,a,b) l_scaledbetapdf(x,a,b,0,1);
  unifdens = @(x,a,b) l_uniformpdf(x,a,b);
  truncLogNormpdf = @(x,mu,sigma,xl,xr) l_truncLogNormpdf(x,mu,sigma,xl,xr);
end

%% load hyperparameters
if isempty(hypfile) % use standard
  c19prior = l_gethyp();
else % use other design
  hyp0 = l_gethyp(); % used to fill empty 
  load(hypfile); % load c19prior
  name0 = fieldnames(hyp0);
  nameC19 = fieldnames(c19prior);
  
  exclude = {'IC_HOSP'};
  missing = name0(~fsetop('ismember',name0,nameC19));
  missing = missing(~fsetop('ismember',missing,exclude));
  for k = 1:numel(missing)
    c19prior.(missing{k}) = hyp0.(missing{k});
  end
end

%% sigma | rate of exposure
if isstruct(X) % probability
  rates.sigma = l_scaledbetapdf(X.sigma,c19prior.sigmaHyp(1),c19prior.sigmaHyp(2),...
    c19prior.sigmaPQ(1), c19prior.sigmaPQ(2));
elseif X < inf % sample
  rates.sigma = c19prior.sigmaPQ(1) + diff(c19prior.sigmaPQ)*...
    betarnd(c19prior.sigmaHyp(1),c19prior.sigmaHyp(2),1,Nreplicas);
else % mean
  rates.sigma =  c19prior.sigmaPQ(1) + diff(c19prior.sigmaPQ)*c19prior.sigmaHyp(1) ...
    / sum(c19prior.sigmaHyp);
end

%% gammaI/A | Symptomatic/Asympomatic waiting time
if fix % fixate dimensions
  rates.gammaI = rates.sigma;
  rates.gammaA = rates.sigma;
elseif isstruct(X) % probability
  rates.gammaI = l_scaledbetapdf(X.gammaI,...
    c19prior.gammaIHyp(1),c19prior.gammaIHyp(2),c19prior.gammaIPQ(1),c19prior.gammaIPQ(2));
  %rates.gammaI = betapdf(X.gammaI, gammaIHyp(1), gammaIHyp(2));
  rates.gammaA = rates.gammaI;
  %rates.gammaA = betapdf(X.gammaA, gammaAHyp(1), gammaAHyp(2));
elseif X < inf % sample
  rates.gammaI = c19prior.gammaIPQ(1) + ...
    diff(c19prior.gammaIPQ)*betarnd(c19prior.gammaIHyp(1),c19prior.gammaIHyp(2),...
    1,Nreplicas);
  %rates.gammaI = betarnd(gammaIHyp(1), gammaIHyp(2), 1, Nreplicas);
  rates.gammaA = rates.gammaI;
  %rates.gammaA = betarnd(gammaAHyp(1), gammaAHyp(2), 1, Nreplicas);
else % mean
  rates.gammaI = c19prior.gammaIPQ(1) + diff(c19prior.gammaIPQ)*c19prior.gammaIHyp(1)/...
    sum(c19prior.gammaIHyp);
  %rates.gammaI = gammaIHyp(1)/sum(gammaIHyp);
  rates.gammaA = rates.gammaI;
  %rates.gammaA = gammaAHyp(1)/sum(gammaAHyp);
end

%% gammaH | hosptial waiting time
if isstruct(X) % probability
  rates.gammaH = l_scaledbetapdf(X.gammaH,c19prior.gammaHHyp(1),...
    c19prior.gammaHHyp(2), c19prior.gammaHPQ(1), c19prior.gammaHPQ(2));
elseif X < inf % sample
  rates.gammaH = c19prior.gammaHPQ(1) + diff(c19prior.gammaHPQ)*...
    betarnd(c19prior.gammaHHyp(1),c19prior.gammaHHyp(2),1,Nreplicas);
else % mean
  rates.gammaH = c19prior.gammaHPQ(1) + diff(c19prior.gammaHPQ)*...
    c19prior.gammaHHyp(1) / sum(c19prior.gammaHHyp);
end

%% gammaW | ICU waiting time
if isstruct(X) % probability
  rates.gammaW = l_scaledbetapdf(X.gammaW, c19prior.gammaWHyp(1),...
    c19prior.gammaWHyp(2), c19prior.gammaWPQ(1), c19prior.gammaWPQ(2));
elseif X < inf % sample
  rates.gammaW = c19prior.gammaWPQ(1) + diff(c19prior.gammaWPQ)*...
    betarnd(c19prior.gammaWHyp(1),c19prior.gammaWHyp(2),1,Nreplicas);
else % mean
  rates.gammaW = c19prior.gammaWPQ(1) + diff(c19prior.gammaWPQ)*...
    c19prior.gammaWHyp(1) / sum(c19prior.gammaWHyp);
end

%% lambda | commuting transport
% commuting transport: scaled beta prior lambda ~ tau * (0.5+B(2,2)),
% with tau = work fraction of time
if isstruct(X) % probability
  rates.lambda = l_scaledbetapdf(X.lambda,c19prior.lambdaHyp(1),...
    c19prior.lambdaHyp(2),c19prior.lambdaScaleShift(1),c19prior.lambdaScaleShift(2));
else
  if X < inf % sample
    rates.lambda = c19prior.lambdaScaleShift(1)+diff(c19prior.lambdaScaleShift)*...
      betarnd(c19prior.lambdaHyp(1),c19prior.lambdaHyp(2),nslab,Nreplicas);
  else % mean
    rates.lambda = repmat(c19prior.lambdaScaleShift(1)+...
      diff(c19prior.lambdaScaleShift)*c19prior.lambdaHyp(1)/sum(c19prior.lambdaHyp),nslab,1);
  end
end

%% half_life/rho |
% survival half-time: uniform prior (1,12) hours
if isstruct(X) % probability
  rates.half_life = unifdens(X.half_life, c19prior.hfPQ(1), c19prior.hfPQ(2));
else
  if X < inf % sample
    rates.half_life = unifrnd(c19prior.hfPQ(1), c19prior.hfPQ(2), 1, Nreplicas);
  else % mean
    rates.half_life = repmat(sum(c19prior.hfPQ)/2,1,1);
  end
  rates.rho = -log(0.5)./rates.half_life;
end

%% thetaA/E
% virus shedding from the different categories (parameter thetaI is
% not free thanks to non-dimensionalization).
if ~isstruct(X) % sample
  rates.thetaI_ = ones(1, Nreplicas);
end

if isstruct(X) % probability
  rates.thetaA_ = l_scaledbetapdf(X.thetaA_, c19prior.thetaAHyp(1),...
    c19prior.thetaAHyp(2),c19prior.thetaAPQ(1),c19prior.thetaAPQ(2));
  rates.thetaE_ = l_scaledbetapdf(X.thetaE_, c19prior.thetaEHyp(1),...
    c19prior.thetaEHyp(2),c19prior.thetaEPQ(1),c19prior.thetaEPQ(2));
else % sample
  if X < inf
    rates.thetaA_ = c19prior.thetaAPQ(1)+diff(c19prior.thetaAPQ)*...
      betarnd(c19prior.thetaAHyp(1),c19prior.thetaAHyp(2),1,Nreplicas);
    rates.thetaE_ = c19prior.thetaEPQ(1)+diff(c19prior.thetaEPQ)*...
      betarnd(c19prior.thetaEHyp(1),c19prior.thetaEHyp(2),1,Nreplicas);
  else % mean
    rates.thetaA_ = c19prior.thetaAPQ(1)+diff(c19prior.thetaAPQ)*...
      c19prior.thetaAHyp(1)/sum(c19prior.thetaAHyp);
    rates.thetaE_ = c19prior.thetaEPQ(1)+diff(c19prior.thetaEPQ)*...
      c19prior.thetaEHyp(1)/sum(c19prior.thetaEHyp);
  end
  rates.thetaI = rates.thetaI_.*rates.rho;
  rates.thetaA = rates.thetaA_.*rates.thetaI;
  rates.thetaE = rates.thetaE_.*rates.thetaI;
end

%% R0
% non-spatial R0
if isstruct(X) % probability
  if nativepdf
    R0Dist = makedist('lognormal',c19prior.R0Hyp(1),c19prior.R0Hyp(2));
    R0Dist = truncate(R0Dist, c19prior.R0PQ(1), c19prior.R0PQ(2));
    rates.R0 = truncLogNormpdf(R0Dist, X.R0);
  else
    rates.R0 = truncLogNormpdf(X.R0,c19prior.R0Hyp(1),c19prior.R0Hyp(2),...
      c19prior.R0PQ(1), c19prior.R0PQ(2));
  end
else
  % only create the distribution if we want a sample
  R0Dist = makedist('lognormal',c19prior.R0Hyp(1),c19prior.R0Hyp(2));
  R0Dist = truncate(R0Dist, c19prior.R0PQ(1), c19prior.R0PQ(2));
  if X < inf % sample
    rates.R0 = random(R0Dist, [nslab,Nreplicas]);
  else % mean
    rates.R0 = repmat(mean(R0Dist),nslab,1);
  end
end

%% IFR
if isstruct(X) % probability
  rates.IFR = l_scaledbetapdf(X.IFR,c19prior.IFRHyp(1),c19prior.IFRHyp(2),...
    c19prior.IFRPQ(1),c19prior.IFRPQ(2));
else % sample
  if X < inf % sample
    rates.IFR = c19prior.IFRPQ(1)+diff(c19prior.IFRPQ)*...
      betarnd(c19prior.IFRHyp(1),c19prior.IFRHyp(2),nslab,Nreplicas);
  else % mean
    rates.IFR = c19prior.IFRPQ(1)+(c19prior.IFRPQ(2)-c19prior.IFRPQ(1))*c19prior.IFRHyp(1)/...
      sum(c19prior.IFRHyp);
    rates.IFR = repmat(rates.IFR,nslab,1);
  end
end

%% E2I | proportion of individuals from E -> I, where 1-F0: E -> A.
if isstruct(X) % probability
  rates.E2I = l_scaledbetapdf(X.E2I,c19prior.E2IHyp(1),c19prior.E2IHyp(2),...
    c19prior.E2IPQ(1), c19prior.E2IPQ(2));
  rates.A2I = ones(size(rates.E2I));
else % sample
  if X < inf
    rates.E2I = c19prior.E2IPQ(1) + diff(c19prior.E2IPQ)*...
      betarnd(c19prior.E2IHyp(1),c19prior.E2IHyp(2),1,Nreplicas);
  else % mean
    rates.E2I = c19prior.E2IPQ(1)+diff(c19prior.E2IPQ)*c19prior.E2IHyp(1)/sum(c19prior.E2IHyp);
  end
  rates.A2I = zeros(size(rates.E2I));
end

%% SIR_MORT
SIR_MORT = c19prior.SIR_MORT;

%% HOSP_MORT
HOSP_MORT = c19prior.HOSP_MORT;

%% IC_HOSP
if isempty(c19prior.IC_HOSPPQ) % legacy version | fixed value
  if isstruct(X) % probability
    rates.IC_HOSP = ones(1,Nreplicas);
  else
    rates.IC_HOSP = repmat(c19prior.IC_HOSPmean,1,Nreplicas);
  end
else % sample IC_HOSP
  if isstruct(X) % probability
    rates.IC_HOSP = l_scaledbetapdf(X.IC_HOSP,c19prior.IC_HOSPHyp(1),...
      c19prior.IC_HOSPHyp(2),...
      c19prior.IC_HOSPPQ(1),c19prior.IC_HOSPPQ(2));
  else % sample
    if X < inf % sample
      rates.IC_HOSP = c19prior.IC_HOSPPQ(1)+diff(c19prior.IC_HOSPPQ)*...
        betarnd(c19prior.IC_HOSPHyp(1),c19prior.IC_HOSPHyp(2),1,Nreplicas);
    else % mean
      rates.IC_HOSP = c19prior.IC_HOSPPQ(1)+(c19prior.IC_HOSPPQ(2)-...
        c19prior.IC_HOSPPQ(1))*c19prior.IC_HOSPHyp(1)/...
        sum(c19prior.IC_HOSPHyp);
    end
  end
end

%% HOSP | upper limit depends on IFR, HOSP_MORT, and IC_HOSP
if isstruct(X) % probability
  rates.HOSP = l_scaledbetapdf(X.HOSP,c19prior.HOSPHyp(1),c19prior.HOSPHyp(2),...
    c19prior.HOSPPQ(1),c19prior.HOSPPQ(2));
else % sample
  rates.HOSP = c19prior.HOSPPQ(2).*betarnd(c19prior.HOSPHyp(1),c19prior.HOSPHyp(2),...
    1,Nreplicas);
  if X < inf % sample
    rates.HOSP = c19prior.HOSPPQ(2).*betarnd(c19prior.HOSPHyp(1),c19prior.HOSPHyp(2),...
      1,Nreplicas);
  else % mean
    rates.HOSP = c19prior.HOSPPQ(2).*c19prior.HOSPHyp(1)/sum(c19prior.HOSPHyp);
  end
end


%% Demographic averages:
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

% beta | deduced contraction rate beta
rates.beta = rates.R0.^2./...
  (rates.thetaE_./rates.sigma+ ...
  (1-rates.F0ave).*rates.thetaA_./rates.gammaA+ ...
  (rates.F0ave + (1-rates.F0ave).*rates.F1ave)./rates.gammaI);

names = fieldnames(rates);
if isstruct(X) % probability of rates
  % compute the sum of the log densities.
  rates.sumlog_dynamic = zeros(nslab,Nreplicas);
  dynamic = {'lambda' 'R0' 'IFR'};
  exclude = {'lambda' 'beta' 'F0ave' 'F1ave' 'F2ave' 'F2dave' 'F3ave' ...
    'F3dave' 'F4ave' 'beta'};
  for j = 1:nslab
    for i = 1:numel(dynamic)
      if ~any(strcmp(dynamic{i},exclude))
        dens = log(rates.(dynamic{i}));
        rates.sumlog_dynamic(j,:) = rates.sumlog_dynamic(j,:) + dens(j,:);
      end
    end
  end
  
  rates.sumlog = 0;
  for i = 1:numel(names)
    if ~any(strcmp(names{i},[exclude dynamic]))
      rates.sumlog = rates.sumlog + log(rates.(names{i}));
    end
  end
  rates.sumlog = rates.sumlog + sum(rates.sumlog_dynamic,1);
end

% .meta
cellrates = struct2cell(rates);
try
  rates.meta.hash = fsetop('check',reshape(cat(1,cellrates{:}),[],1));
catch
  1;
end
rates.meta.rev = c19prior.rev;
% ----------------------------------------------------------------------
function dens = l_scaledbetapdf(x, A, B, p, q)
% DENS = L_SCALEDBETAPDF density of the scaled beta distribution.
%
% DENS = L_SCALEDBETAPDF(X,A,B,P,Q) according to which:
%   X = value to evaluate.
%   A = the beta "shape" parameter
%   B = the beta "scale" parameter, see betapdf
%   P = the lower (left) limit of the scaled density
%   Q = the upper (right) limit of the scaled density

nomin = (x-p).^(A-1) .* (q - x).^(B-1);
denom = (q - p).^(A+B-1) .* beta(A, B);
dens = nomin./denom;

dens(x < p) = 0;
dens(x > q) = 0;

% ----------------------------------------------------------------------
function dens = l_truncLogNormpdf(x,mu,sigma,xl,xr)
%DENS = L_TRUNCLOGNORM(X,MU,SIGMA,XL,XR) computes the pdf of a truncated
%   lognormal distribution, with the aim of being faster than makedist
%   route.
%
%   X - point to evalutate
%   MU - normal mean
%   SIGMA - normal standard deviations
%   XL - left truncation
%   XR - right truncation
%
%   Implementation from
%   http://www.isaacpub.org/images/PaperPDF/AdAp_100070_2017041816215848323.pdf

mu = exp(mu);
a = 0.5*sqrt(2)/sigma*log(xl/mu);
b = 0.5*sqrt(2)/sigma*log(xr/mu);

dens = (sqrt(2)*exp(-0.5*(1/sigma^2)*(log(x/mu)).^2)) ./ ...
  (-sqrt(pi)*sigma*(erf(a) - erf(b))*x);

dens(x < xl) = 0;
dens(x > xr) = 0;

% ----------------------------------------------------------------------
function dens = l_uniformpdf(x,xl,xr)
%DENS = L_UNIFORMPDF(X,XL,XR) computes the pdf of the uniform distribution.
%
%   X - point to evaluate
%   XL - left border
%   XR - right border

dens = zeros(size(x));
inside = (xl <= x) & (x <= xr);
dens(inside) = 1/(xr-xl);

% ----------------------------------------------------------------------
function hyp0 = l_gethyp()
%HYP0 = L_GETHYP return the standard hyperparameters in a struct form.

%% sigma
% https://www-acpjournals-org.ezproxy.its.uu.se/doi/pdf/10.7326%2FM20-0504
% suggests 1/sigma ~ Gamma(shape=5.807, scale=0.948);
sigmaPQ = [0.1 0.5]; % prior support
sigmamean = 1/5; % prior mean
sigmaHyp(1) = 2;
sigmaHyp(2) = sigmaHyp(1)/((sigmamean-sigmaPQ(1))/diff(sigmaPQ))-sigmaHyp(1);

%% gammaI/gammaA
gammaIPQ = [0.1 0.5]; % prior support
gammaImean = 1/7; % prior mean

gammaIHyp(1) = 2;
gammaIHyp(2) = gammaIHyp(1)/((gammaImean-gammaIPQ(1))/diff(gammaIPQ))-gammaIHyp(1);

%% gammaH
% See statistics about hospital waiting times in /data/sources/SS/statistik-covid19-vardforlopp.xlsx
% Computed empirical mean is 1/8.26 days^-1
% For computation, see /data/sources/findGamma.R
gammaHPQ = [1/9 1/6]; % prior support
gammaHmean = 1/8.26; % prior mean 0.1211

gammaHHyp(1) = 2;
gammaHHyp(2) = gammaHHyp(1)/((gammaHmean-gammaHPQ(1))/diff(gammaHPQ))-gammaHHyp(1);

%% gammaW
% See statistics about ICU waiting times in /data/sources/SS/statistik-covid19-vardforlopp.xlsx
% Computed empirical mean is 1/11.89 days^-1
% For computation, see /data/sources/findGamma.R
gammaWPQ = [1/15 1/10];
gammaWmean = 1/11.89; % prior mean 0.0841;%

gammaWHyp(1) = 2;
gammaWHyp(2) = gammaWHyp(1)/((gammaWmean-gammaWPQ(1))/diff(gammaWPQ))-gammaWHyp(1);

%% lambda
tau = 5/7*8/24;
lambdaHyp = [2, 2];
lambdaPQ = tau*[0.5, 1.5];

%% half_life
hfPQ = [1 12]/24;

%% thetaA/E
thetaAHyp = [2 2];
thetaEHyp = thetaAHyp;
thetaAPQ = [0 2];
thetaEPQ = thetaAPQ;

%% R0
R0Hyp = [log(1.3) 0.4];
R0PQ = [0 4];

%% IFR
IFRHyp = [2 4];
IFRPQ = [0 0.02]; % supp(ifr) = [0,2]%
IFRmean = [];

%% F0
E2IPQ = [0 1]; % prior support
E2IHyp(1) = 44.1089;
E2IHyp(2) = 15.6984;

%% SIR MORT
SIR_MORT = 0.2031;

%% HOSP MORT
HOSP_MORT = 0.4918;

%% IC_HOSP
IC_HOSPPQ = [];
IC_HOSPHyp = [];
IC_HOSPmean = 0.1151;

%% I_MORT
I_MORTmean = [];

%% HOSP
% compute the I_IFR_scaling, and mean for HOSP.
% we need this for IFR below:
E2Imean = E2IPQ(1) + diff(E2IPQ)*E2IHyp(1) / sum(E2IHyp);
A2Imean = 0;

I_IFR_scale = 1/(E2Imean + (1 - E2Imean)*A2Imean);
I_IFR_mean = (IFRPQ(1)+(IFRPQ(2)-IFRPQ(1))*IFRHyp(1)/sum(IFRHyp)) * I_IFR_scale;

HOSPHyp = [2 2];
HOSPmean = [];
HOSPPQ = [0 I_IFR_mean./SIR_MORT ... % upper boundary
  ./(HOSP_MORT+IC_HOSPmean) ...
  .*(1-IC_HOSPmean*(1-SIR_MORT))];

%% pack into struct
hyp0 = struct();
% pack up as an output struct
hyp0.sigmaPQ = sigmaPQ;
hyp0.sigmaHyp = sigmaHyp;
hyp0.sigmamean = sigmamean;

hyp0.gammaIPQ = gammaIPQ;
hyp0.gammaIHyp = gammaIHyp;
hyp0.gammaImean = gammaImean;

hyp0.gammaHPQ = gammaHPQ;
hyp0.gammaHHyp = gammaHHyp;
hyp0.gammaHmean = gammaHmean;

hyp0.gammaWPQ = gammaWPQ;
hyp0.gammaWHyp = gammaWHyp;
hyp0.gammaWmean = gammaWmean;

hyp0.lambdaPQ = lambdaPQ;
hyp0.lambdaHyp = lambdaHyp;
hyp0.lambdaScaleShift = lambdaPQ;

hyp0.hfPQ = hfPQ;

hyp0.thetaAPQ = thetaAPQ;
hyp0.thetaAHyp = thetaAHyp;

hyp0.thetaEPQ = thetaEPQ;
hyp0.thetaEHyp = thetaEHyp;

hyp0.R0PQ = R0PQ;
hyp0.R0Hyp= R0Hyp;

hyp0.E2IPQ = E2IPQ;
hyp0.E2IHyp = E2IHyp;
hyp0.E2Imean = E2Imean;

hyp0.IFRPQ = IFRPQ;
hyp0.IFRHyp = IFRHyp;
hyp0.IFRmean = IFRmean;

% sole two point estimates:
hyp0.SIR_MORT = SIR_MORT;
hyp0.HOSP_MORT = HOSP_MORT;

hyp0.IC_HOSPPQ = IC_HOSPPQ;
hyp0.IC_HOSPHyp = IC_HOSPHyp;
hyp0.IC_HOSPmean = IC_HOSPmean;

hyp0.HOSPPQ = HOSPPQ;
hyp0.HOSPHyp = HOSPHyp;
hyp0.HOSPmean = HOSPmean;

% sampled later on the fly, but for convenience:
hyp0.I_MORTmean = I_MORTmean;


hyp0.rev = '26-May-2021'; % latest day the prior was updated

% ----------------------------------------------------------------------