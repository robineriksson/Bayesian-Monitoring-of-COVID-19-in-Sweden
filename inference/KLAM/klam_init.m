function [thetas,sl,slab,amparam,amfunc,outverb,Ydata] = klam_init(...
  region,nMCMC,verb,nslab,date,init,scaleS,register,useCSSS, ...
  perslab,fix,ksmc,state0)
%KLAM_INIT initializes the Kalman Likelihood Adaptive Metropolis
% sampler. The object is constructed with the data (YDATA) and so on. Finally,
% the AM sampler is called and NMCMC samples are generated from the posterior.
%
% KLAM_INIT computes the (approximate) Bayesian posterior of the specified
%   REGION until the date DATE, using NMCMC samples in an Adaptive
%   Metropolis (AM) chain. Init holds initial guess variables for the AM
%   chain, however, this can be left empty for default values. OUTVERB
%   gives diagnostic information from the AM solver, e.g., acceptance_rate,
%   prior_failed, earlyReject, stopper_slabs, tt.
%
%   REGION - the data from which region to use in the inference, e.g.,
%   'Stockholm' or 'Uppsala'. If all data use  '', if aggregated Sweden
%   without transport use 'Sweden'.
%
%   NMCMC - The length of the Adaptive Metropolis chain. Consider that you
%   should include a burn-in distance, which can be found post the
%   generation of the chain, i.e., has the chain entered a 'stable' region?
%
%   VERB - specifies if debugging messages should be displayed
%
%   NSLAB - defines the number of slabs (time periods) to use. If INF, they
%   are automaticaly generated as per month, with an expection to the final
%   slab which will be the penultimate month + 13 days into the ultimate
%   month if the final date is < 14th in the month, as a minimum of 14 days
%   per slab is set.
%
%   DATE - [startdate enddate] the start and stop of loaded data.  If
%   DATE(1) < Date.date(1) then DATE(1) <- Data.date(1). And if
%   DATE(2) > Date.date(end) then DATE(2) <- Data.date(end).
%
%   INIT - empty or struct with fieldnames *rates0* which is the inital
%   guess for the chain, *ratesAll* which is all the rates in the posterior
%   to be used for the covariance estimate for the AM chain, and *nslab*
%   number of slabs used for the initial guess.
%
%   SCALES - step scaling parameter Sd affects the acceptance rate, in
%   general, smaller Sd gives higher acceptance rate and vice versa. SCALES
%   then allows to scales the standard Sd implemented in this function,
%   also visible in the struct amparam which holds method parameters.
%
%   REGISTER - what main data source (can only be one), e.g., 'C19',
%   'C19pred', 'FHM'.
%
%   USECSSS - add an estimate of the I compartment.
%
%   PERSLAB - if true then use AM_SEC.m else AM.m is used.
%
%   FIX     - if sigma = gammaI/A or not
%
%   KSMC     - true if SMC sampler, else AM sampler
%
%   STATE0 - initial state to start simulations from. If left empty,
%   start from scratch.
%
%   See also, AM.m which is called from this function to generate the
%   posterior samples.

% R. Eriksson 2020-10-06

switch nargin
  case 1
    nMCMC=1e4;
    verb=false;
    nslab=inf;
    date=[200319, 201231];
    init=[];
    scaleS=1;
    register='C19';
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 2
    verb=false;
    nslab=inf;
    date=[200319, 201231];
    init=[];
    scaleS=1;
    register='C19';
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 3
    nslab=inf;
    date=[200319, 201231];
    init=[];
    scaleS=1;
    register='C19';
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 4
    date=[200319, 201231];
    init=[];
    scaleS=1;
    register='C19';
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 5
    init=[];
    scaleS=1;
    register='C19';
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 6
    scaleS=1;
    register='C19';
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 7
    register='C19';
    useCSSS=false;
    fix=false;
    ksmc=0;
    state0 = [];
  case 8
    useCSSS=false;
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 9
    perslab=true;
    fix=false;
    ksmc=0;
    state0 = [];
  case 10
    fix=false;
    ksmc=0;
    state0 = [];
  case 11
    ksmc=0;
    state0 = [];
  case 12
    state0 = [];
  case 13
    ; %do nothing
  otherwise
    error('wrong number of input parameters');
end

if isempty(state0)
  state0 = struct('Y0', [], 'Y0Cov', [], 'date', [], 'file', []);
end
%% Observation model
%
% The parameterization is performed on a single region, however as the code
% is constructed for a network of connected regions, we need to account for
% that and supply the correct parameters.

% specify observation model
obsrates.states = {5, 6, 7}; % state #5, #6, #7
obsrates.indobs = {}; % No indirect measurements
obsrates.indobspars = {};
obsrates.nstate = 8;
obsrates.R0 = 1*ones(3,1);
obsrates.rdiag = (1e-3)^2*ones(3,1);


% transmission matrix D
load Dcounties % D, lan


%% Prepare data
% load data
if contains(register,'URDME') % synthetic data
  filename = mfilename('fullpath');
  Data_raw = load([filename(1:end-24) ...
    'URDME/URDMEoutput/URDME_all']);

  Data = struct();
  regid = strcmp(Data_raw.D.regions,region);
  Data.regions = Data_raw.D.regions(regid);
  Data.date = Data_raw.D.date;
  % which data "page"
  try
    runid = str2double(register(end));
  catch
    runid = 1;
  end
  Data.H = Data_raw.D.U(:,5,regid,runid);
  Data.W = Data_raw.D.U(:,6,regid,runid);
  Data.D = Data_raw.D.U(:,7,regid,runid);

  Data.hash = fsetop('check',Data_raw.D.U(:));

  Data.rev = Data_raw.D.date(end);
  Data.reg = register;
  T = 1;
  D = 1;
  lan = {Data.regions};

else % load empirical data
  Data = loadData(register);
  Data = polishData(Data,'D','Dinc',1);
  Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

  if strcmp(region, 'Sweden') % Sweden is aggregated, no transport
    lan = {'Sweden'};
    Data.H = sum(Data.H,2);
    Data.W = sum(Data.W,2);
    Data.D = sum(Data.D,2);
    T = 1;
    D = 1;
  elseif any(strcmp(region, Data.regions))% other: Stockholm, Uppsala, etc.
    ixRegion = find(strcmp(region,Data.regions));
    lan = lan(ixRegion);
    fnames = fieldnames(Data);
    for k = 1:numel(fnames)
      name = fnames{k};
      if ~any(strcmp({'date', 'hash', 'rev', 'reg'},name))
        Data.(name) = Data.(name)(:,ixRegion);
      end
    end
    T = 1;
    D = 1;
  elseif strcmp(region, 'super1')
    regions = {'Uppsala' 'Dalarna' 'Gävleborg' 'Västmanland'};

    ixRegion = find(fsetop('ismember',Data.regions,regions));

    lan = {reshape(char(regions)',1,[])};
    fnames = fieldnames(Data);
    for k = 1:numel(fnames)
      name = fnames{k};
      if ~any(strcmp({'date', 'hash', 'rev', 'reg'},name))
        Data.(name) = Data.(name)(:,ixRegion);
        if ~strcmp(name,'regions') && size(Data.(name),2) > 1
          Data.(name) = sum(Data.(name),2);
        end
      end
    end
    T = 1;
    D = 1;

  else % use all data, no aggregate, i.e., Sweden with transport
    T = [speye(numel(lan)); sparse(1,1:numel(lan),1)];
  end

end



% keep information about the data in the posterior
amparam.dataHash = Data.hash;
amparam.dataRev = Data.rev;
amparam.dataReg = Data.reg;


% upper and lower bound on the dates.
% If the dates supplied aren't included in the data, return the maxrange.
if date(1) > date(2)
  error('Date(1) > Date(2) | not allowed');
end
if date(1) < Data.date(1)
  warning("Supplied starting date < Data.date(end), overwritten")
  date(1) = Data.date(1);
end
if date(2) > Data.date(end)
  warning("Supplied stopping date > Data.date(end), overwritten")
  date(2) = Data.date(end);
end

if ~isempty(state0.file)
  date(1) = state0.date(1);
  amparam.origo = state0; % pointer reference for future merge.
end


% concatenate the data into an object which is easier to handle.
ixdata = find(Data.date == date(1)):find(Data.date == date(2));
Ydata = cat(3,Data.H(ixdata,:),Data.W(ixdata,:),Data.D(ixdata,:));
Ydata = permute(Ydata,[3 2 1]);

amparam.date = Data.date(ixdata); % save dates

% Note, CSSS data is not included in the release as it was not used in
% the the paper to achieve any of the presented result.
if useCSSS
  % define a common time frame, including dates
  TSPAN = Data.date;
  DATES = datenum(2e3+floor(Data.date(1)/1e4), ...
    mod(floor(Data.date(1)/1e2),1e2) ,...
    mod(Data.date(1),1e2));
  DATES = DATES:DATES+numel(TSPAN)-1;
  DATES = datevec(DATES);
  DATES = DATES(:,1:3)*[1e4 1e2 1]'-2e7; % finally!


  % CSSS-data
  csss = loadData('CSSS');
  Ihi = csss.Ihi(:,ixRegion);
  Ilo = csss.Ilo(:,ixRegion);
  Imid = csss.Imid(:,ixRegion);

  Ihi(csss.date > 201231) = []; % data post 31 Dec 2020, is "bad".
  Imid(csss.date > 201231) = [];
  Ilo(csss.date > 201231) = [];
  csss.date(csss.date > 201231) = [];

  ixcsss = find(DATES == csss.date(1)):find(DATES == csss.date(end));

  % settle on a relative variance from given asymmetric 95% CI
  % Assume a larger error than what they present.
  varI = mean(((Ihi-Ilo)/2./Imid))^2;

  % insert CSSS data
  Idata = permute([NaN(ixcsss(1)-ixdata(1),1); Imid; NaN(ixdata(end)-ixcsss(end),1)],[3 2 1]);
  Ydata = cat(1,Ydata,Idata);
  %
  amparam.dataHash = [amparam.dataHash csss.hash];
  amparam.dataRev  = [amparam.dataRev ' | ' csss.rev];
  amparam.dataReg  = [amparam.dataReg ' | ' csss.reg];


  % demographics
  load Ncounties
  Npop = sum(N,1);

  % correction factor in xmod = fac*x
  load Fcases_interp
  fac = sum(Eave(3:end))/sum(Iave(3:end));

  % specify observation model
  obsrates.indobs = {1}; % CSSS data for the I-compartment
  obsrates.indobspars = {100/(Npop(ixRegion)*fac)};
  obsrates.R0 = [ones(3,1); 1e-4];
  obsrates.rdiag = [(1e-3)^2*ones(3,1); varI];
end



%% specify output model
H = getC19obs(obsrates); % i.e., "H", "W", "D", ["I (in percent)" if csss]
% read as: "Stockholm", "Uppsala", "Sweden total"
G = kron(T,H);

% uncertainty parameters
Q.Q0 = speye(numel(lan)*obsrates.nstate);
Q.qdiag = 0.05^2;



%% Slabs
% define the slabs
amparam.nslab=nslab;
if ~isinf(nslab)
  % if not, we generate the stop dynamically.
  slabstop = zeros(amparam.nslab+1,1);
end

%
switch amparam.nslab
  case 1 % [1,1) per month
    slabstart = (mod(Data.date(ixdata),100) == 1)';
  case 8 % [8,8) per month
    slabstart = (mod(Data.date(ixdata),100) == 8)';
  case 15 % [15, 15) per month
    slabstart = (mod(Data.date(ixdata),100) == 15)';
  case 22 % [22, 22) per month
    slabstart = (mod(Data.date(ixdata),100) == 22)';
  otherwise
    error("the slab case is not defined. Try any of [1 8 15 22]")
end
slabstart_find = find(slabstart);

% first and last period might be small, remove if they are
if slabstart_find(1) < 13
  disp('first period too small, merge with second');
  slabstart(slabstart_find(1)) = 0;
end
if numel(slabstart) - slabstart_find(end) < 13
  disp('final period too small, merge with second to last');
  slabstart(slabstart_find(end)) = 0;
end
slabstop = [1 find(slabstart) numel(ixdata)];
[~,~, slab] = histcounts(Data.date(ixdata), [Data.date(ixdata(slabstart)); inf]);
slab=(slab+1)'; % zero indexed.

if numel(slab) ~= numel(ixdata)
  error("slab and data is not the same length")
end


amparam.nslab = numel(slabstop)-1;
amparam.slabstop = slabstop;


%% Initial guess and ratenames
%
% amparam  which was introduced earlier, holds AM specific parameters
% that are needed. These range from input to the likelihood function to
% the sample function.
%
% Modelspecific parameter, see priorenger for further usage.
amparam.fix   = fix; % sigma = gammaI = gammaA

% AM-debug parameters
savepath = mfilename('fullpath');
amparam.savepath = [savepath(1:end-14) 'results/' region]; % where to save intermediate results.
if isempty(verb)
  verb = true;
end
amparam.verb = verb; % display debug information or not in AM.


% Names of the rates that are to be inverted using the AM. The names
%  cover the fixed and the dynamic (slab dependent ones)
%
% The 'fixed' rates are the ones that do not change over slabs.
amparam.fixed = {'sigma' 'gammaH' 'gammaW' ...
  'thetaA_' 'thetaE_' ...
  'E2I' 'half_life' 'HOSP' 'IC_HOSP'};
if ~amparam.fix
  amparam.fixed = [amparam.fixed 'gammaI'];% 'gammaA'];
end
%  While the dynamic ones are changing values over slabs.
amparam.dynamic = {'R0' 'IFR'};
if numel(D) > 1
  amparam.dynamic = [amparam.dynamic 'lambda'];
end
amparam.dynamic = repmat(amparam.dynamic, amparam.nslab, 1);
amparam.dynamic = amparam.dynamic(:)';


% concatenate all the fieldnames.
amparam.ratenames = [amparam.fixed amparam.dynamic];

% get hyperparameters for prior.
amparam.hypfile = [savepath(1:end-24) 'inference/c19prior'];
[~,~,amparam.hyp] = priorenger(inf,amparam.fix,amparam.nslab,...
  amparam.hypfile);


% initial guess
if ksmc % use SMC
  rates = priorenger(nMCMC, amparam.fix, amparam.nslab,amparam.hypfile);
else % use AM
  if isstruct(init)

    r0 = init.rates0; % extract prior initial knowledge
    r0names = init.names; % make sure that we unpack the correct names.
    ratesAll = init.ratesAll; % will need below.

    % if we have an inital state, then we should remove the dynamic ones,
    % that aren't needed.
    if ~isempty(state0.file)
      % this case, we have an inital state file, and a inital posterior
      % file. It could be, that the inital posterior file was recorded for
      % a longer time period than the one that we want now. For this, we
      % then have to remove and only use the parameters that covers the
      % periods that we are of interest now, e.g., posterior file is
      % generated on data from 200401 to 210501, and we want to produce a
      % posterior from 210401 to 210601. In this scenario we want to
      % exclude dynamic slabs that are before 210401 but the ones after
      % could be used as inital guesses.


      % how many slabs to keep
      start = find(init.meta.date  == date(1));
      if false %init.meta.slabstop(end) == start
        keepnum = 0;
      else
        keepnum = numel(find(init.meta.slabstop >= start));
        while keepnum > amparam.nslab % to many guesses.
          keepnum = keepnum + 1;
        end
      end

      % which ones has size which changes dependent on the the number of slabs?
      names_unique = unique(r0names);
      slabdependent = {};
      for i = 1:numel(names_unique)
        name = names_unique{i};
        if numel(find(strcmp(name,r0names))) > 1
          slabdependent{numel(slabdependent)+1} = name;
        end
      end

      removeindex = [];
      for name = slabdependent
        ismem = find(fsetop('ismember',r0names,name));
        removeindex = [removeindex, ismem(1:end-keepnum)];
      end

      amparam.initslabs =init.meta.date(init.meta.slabstop(1:end-keepnum));
      amparam.initslabs(end) = max(amparam.initslabs(end),date(1));

      r0(removeindex) = [];
      r0names(removeindex) = [];
      ratesAll(removeindex,:) = [];
    else
      keepnum = init.nslab;
      start = 0; % not needed.
    end % else keep all previous values.
    rates = mat2struct(r0, r0names, amparam.hyp, amparam.fix, keepnum);


    % extend inital guess if necessary there is a possibility that our
    % prior guess isn't incStockholmluding all the necessary slabs that we want,
    % say for example that this new week that we added a new slab, then
    % the prior initial guess needs to be extended.
    extendFlag=max(0,amparam.nslab - size(rates.R0,1));
    if size(rates.R0,1) < amparam.nslab % extend!
      % recognizes that we need to extend the initial covariance guess as well.
      % generate a sample at the mean of prior
      rates_ = priorenger(inf, amparam.fix, amparam.nslab,amparam.hypfile);

      % merge newly generated sample and the provided sample s.t., the new
      % samples are appended to the end of the provided samples.
      %
      % Example, structMerge does the following for all the ratenames.
      % rates.R0 = [rates.R0; rates_.R0]

      if init.meta.slabstop(end) == start % only store non-dynamic ones
        rates = structMerge(rates,rates_, ...
          {0,1}, {1:size(rates_.R0,1),1});
      else
        rates = structMerge(rates,rates_, ...
          {1:size(rates.R0,1),1}, {keepnum+1:amparam.nslab,1});
      end


    elseif size(rates.R0,1) > amparam.nslab % reduce
      error(['You`ve submitted a initialization file which' ...
        ' hold more slabs than the data you`re solving for has' ...
        ', i.e., you are using future data. NOT allowed.']);
    end

    rates.meta = init.meta;
  else % no prior initialization
    rates = priorenger(inf, amparam.fix, amparam.nslab,amparam.hypfile);
  end
end


% The inference using AM is done in matrix form, and not struct for
% that C19filt is using. Therefore,
% 1. convert struct to a matrix.
[theta0, theta0_names_, meta] = struct2mat(rates);

% 2. in the matrix, remove the rates that we do not want to infer by
% AM, i.e., amparam.ratenames.
theta0 = theta0(fsetop('ismember',theta0_names_,amparam.ratenames),:)';

% 3. adjust the order of the names in amparam.ratenames, s.t. they
% are alligned with the order of the names in the matrix form.
theta0_names = theta0_names_(fsetop('ismember',theta0_names_,amparam.ratenames));
amparam.ratenames = theta0_names;

% posterior samples are stored in the matrix format .
if ksmc
  thetas = zeros([size(theta0),ksmc]);
  thetas(:,:,1) = theta0;
else % AM
  thetas = zeros(numel(theta0),nMCMC);
  thetas(:,1) = theta0;
end


%% log-likelihood function
% The log-likelihood function, for ease of use, should only take input
% that will change durring the iterations. Therefore, we create an
% anonymous function which reduce the number of input needed per
% iteration.

% these will be used durring inference
exception.LB = -1e2;
exception.UB = 1e7;
exception.SDFAC = 0.1;
exception.AbsMagn = 1e3;
amparam.exception = exception;

% and these for the initial sample (to ensure that we can start from
% random)
exception0.LB = -inf;
exception0.UB = inf;
exception0.SDFAC = inf;
exception0.AbsMagn = inf;




if perslab
  logKL_func = @(x,state0,ixslab) logKL(x,state0,G,D,obsrates,Q,Ydata,...
    slab,ixslab,exception);
  [~,~,L] = logKL(rates,state0,G,D,obsrates,Q,Ydata,...
    slab,ixslab,exception0);
else
  logKL_func = @(x,Ydata)logKL(x,state0,G,D,obsrates,Q,Ydata,...
    slab,1:slabstop(end),amparam.exception);
  [~,~,L] = logKL(rates,state0,G,D,obsrates,Q,Ydata,...
    slab,1:slabstop(end),exception0);
end

if ~ksmc && any(isinf(sum(L,1)))
  disp(sum(L)');
  error('initial guess is invalid, Inf likelihood');
end

%% likelihood storage
% We want to keep track of the likelihood of all the samples in the
% Markov chain. Not only for evaluation, but it also alows us to
% continue the chain after stopping.
if ksmc
  sl = sum(L,1);
else
  sl = zeros(size(L,1), nMCMC);
  sl(:,1) = sum(L,2);
end

%% Adaptivity - mean
% The adaptivity of AM requires the mean of the posterior samples to
% be generated at each step, which is done recursively. To initalize
% the adaptivity there is here two alternative 1) a supplied initial
% prior guess -> extract the mean of the samples or 2) the rolling
% mean starts with the inital value in the chain.

if isstruct(init)

  % given all the posterior samples, it can also include rates that are
  % not to be determined by AM. Therefore, we make sure to select only
  % the ones that are given in amparam.ratenames.
  old = ratesAll(fsetop('ismember',r0names,amparam.ratenames),:);
  if(extendFlag)
    % find the indices where the new slabs should be added (exclude these)
    idXtend = zeros(numel(amparam.ratenames),1);
    for name = unique(amparam.dynamic)
      ind = find(strcmp(amparam.ratenames,name));
      idXtend(ind(end-extendFlag+1:end)) = 1;
      %find(strcmp(amparam.ratenames,name),1,'last')) = 1;
    end
    xbar = theta0;
    xbar(~idXtend) = mean(old,2)';
  else
    xbar = mean(old,2)';
  end
  amparam.xbar = xbar;
else % no init
  amparam.xbar = theta0;
end

%% Scaling of AM
% The proposal in AM are scaled by a parameter Sd, and also by a
% paramter e which is there to safeguard from singularity.
amparam.Sd = 2.4^2/numel(theta0); % optimal from Gelman et al (1996)
% increase/decrease the acceptance rate (found from tuning)
amparam.Sd = scaleS*amparam.Sd;
amparam.e = amparam.Sd/100; % e << S

%% Adaptivity - covariance
% AM will start learning the sample covariance after an inital
% burn-in period of i0. Larger i0 signals higher trust in the initial
% covariance
if isstruct(init)
  amparam.i0 = 1e4; % high trust
else
  amparam.i0 = 1e1; % low trust
end

% if supplied initial guess (prior posterior)
if ~ksmc
  if isstruct(init)

    % if need to extend by additional slabs.
    if(extendFlag)
      % create a new eye matrix of desired size.
      ccov = eye(numel(idXtend));
      % add the old covariance matrix where the old slabs are know
      ccov(~idXtend,~idXtend) = cov(old');
      amparam.xcov = ccov;
    else
      % find the covariance of the supplied samples
      amparam.xcov = cov(old');
    end
    % scale is, as would be done in AM.
    amparam.sigma = amparam.Sd*amparam.xcov + amparam.Sd*amparam.e*eye(numel(theta0));
  else % no init

    % empty initalization can be done by a scaled eye matrix or a random one.
    % What one needs to remember is that the covariance matrix is SPD.
    amparam.sigma = 1e-3*eye(numel(theta0));
    amparam.xcov = amparam.sigma;
  end
end





%% Prior
% The acceptance probability needs to consider the prior density of
% the sample, otherwise it would not be truly Bayesian inference.
%
% In our C19-priorenger model, we consider some fixed interval rates,
% meaning, that their density is zero outside this interval. The prior
% density then becomes an important tool to determine when samples
% should be discared before "tested".
%
% As for the log-likelihood function, we supply a reduced
% anom. function for the prior.
prior_func = @(x) priorenger(x,amparam.fix,amparam.nslab,amparam.hypfile);
amparam.prior = prior_func(rates); % *** same here, multiple proposals?
if any(isinf(getfield(amparam.prior,'sumlog')))
  error('initial proposal is rejected based on prior');
end

% get prior revision.
amparam.priorRev = amparam.prior.meta.rev;


% mat2struct in a more compact fashion
m2s = @(x) mat2struct(x, amparam.ratenames, amparam.hyp,...
                      amparam.fix, amparam.nslab);
% the reduced functions are stored in a struct for easy retrieval.
amfunc = struct();
amfunc.logKL = logKL_func;
amfunc.prior = prior_func;
amfunc.m2s = m2s;

%% Generatesamples
if ksmc % use smc sampler
  amfunc.metro = @(nMCMC,thetas,sl,amparam,amfunc,ksmc,empty) ...
    SMC(nMCMC,thetas,sl,amparam,amfunc,ksmc,empty);
else
  if perslab
    amfunc.metro = @(nMCMC,thetas,sl,amparam,amfunc,start,stop) ...
      AM_seq(nMCMC,thetas,sl,amparam,amfunc,start,stop);
  else
    amfunc.metro = @AM;
  end
  ksmc=2; % flag that indicates the solve number of inputs.
end
tic;
[thetas, sl, amparam, outverb] = amfunc.metro(nMCMC, thetas, sl, ...
  amparam, amfunc,Ydata,ksmc,nMCMC);
outverb.tt = toc;
amparam.region = region;
end