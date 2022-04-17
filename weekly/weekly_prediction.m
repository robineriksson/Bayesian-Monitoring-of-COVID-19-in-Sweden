%WEEKLY_PREDICTION Script for generating predictions.
%
%   Either run by copy-and-past from snippets of in the beginning
%   of the file, or run as a script by editing some options:
%
%   test = integer, definest the region
%     0: Sweden
%     4: Stockholm
%     5: Uppsala
%     <many more>
%
%   type =
%     1: without historic lag, but prediction lag days ahead
%     2: with historic lag
%
%   register = string, source of data
%
%   savetofile = {0 1} save to .tex-files
%
%   useCSSS = {0 1} use CSSS-data (experimental)
%
%   phiplotting/Rplotting = {0 1} plot phi and R(ecovered)

% S. Engblom 2021-04-01
% S. Engblom 2020-10-25 (Revision, publication plots)
% S. Engblom 2020-09-18 (Major revision)
% S. Engblom 2020-09-07 (Revision, spatially connected model)
% S. Engblom 2020-09-02 (Revision)
% S. Engblom 2020-04-13

if 0
  % snippets of code for copy and paste (for weekly report)
  % defaults
  posteriordate = '210531'; % change here
  ending        = '1_100'; % at what date does the slabs start.
  datadate     = [200401 210531] % change end date here

  register   = 'C19'; % Employs C19 if not Uppsala, then Uppsala Region
  useCSSS      = false;
  saveall      = true;
  % run all regions
  type = 1, lag = 14;
  regionList = regions(false);

  for reg = [1:21]
    region = regionList{reg}
    FINALRUN = true;
    weekly_prediction
  end

  % lag7-plot, Uppsala
  FINALRUN = true;
  type = 2, lag = 7;
  region = 'Uppsala'
  weekly_prediction

  % lag7-plot, Stockholm
  FINALRUN = true;
  type = 2, lag = 7;
  region = 'Stockholm'
  weekly_prediction

  % lag14-plot, Sweden
  FINALRUN = true;
  type = 1, lag = 14;
  region = 'Sweden'
  weekly_prediction

  % lag7-plot, Sweden
  FINALRUN = true;
  type = 2, lag = 7;
  region = 'Sweden'
  weekly_prediction

  disp('*** Done ***');
end
% *** change here BEGIN { **********************************************
if ~exist('FINALRUN','var') || ~FINALRUN
  % type =
  %   1: without historic lag, but prediction lag days ahead
  %   2: with historic lag
  type = 1
  if type == 1
    lag = 14;
  else
    lag = 7;
  end

  % use CSSS?
  useCSSS=false;

  % save all states
  saveall=true;

  % for reference
  region = 'Uppsala'

  % posterior file date
  posteriordate = '210531';
  ending         = '1_100'; % at what date does the slabs start.

  % dates that "date" the data used
  datadate     = [200401 210531]; % change end date here

  % exclude the network (only affects 'Sweden')
  noNetwork = 1;

  register = 'C19';
else
  disp('*** Running as script... ***');
  FINALRUN = false;
  % (will only run once as script to avoid misunderstandings)
end

% folder for posteriors
abspath = mfilename('fullpath');
prefix = [abspath(1:end-24) 'inference/results/KLAM/'];

% Population data
load Ncounties
Npop = sum(N,1);


% convert region into "test" for data extraction, later.
regionList = regions(false);
test = find(strcmp(region,regionList)) + 3;
if isempty(test)
  if strcmp(region,'Sweden')
    test = 0;
    % construct the Sweden posteriorfile
    if ~exist([prefix 'slam' posteriordate '_Sweden_monthly_' ending '.mat'],'file')
      regs = 1:21;
      files = cell(size(regs));
      for r = 1:numel(regs)
        files{r} = [prefix 'perRegion/slam' posteriordate '_' regionList{regs(r)} '_monthly_' ending '.mat'];
      end
      Weights = Npop;
      rates = posteriorenger(100,files,Weights);
      save([prefix 'slam' posteriordate '_Sweden_monthly_' ending '.mat'],'rates')
    end
  else
    error(['Missing region: ' region]);
  end
end


if test > 0
  posterior = 'perRegion/';
else
  posterior = '';
end

%posterior = [posterior '/smc' posteriordate '_22_' region '_monthly'];
%posterior = [posterior '/slam' posteriordate '_22_' region '_monthly'];
posterior = [posterior 'slam' posteriordate '_' region '_monthly'];

if useCSSS
  posterior = [posterior '_csss'];
end

posterior = [posterior '_' ending '.mat'];

% sample rates
rng(0);
rates = posteriorenger([],[prefix posterior]);
try
  slabs = rates.meta.slab';
catch
  slabs = rates.meta.slabs;
end
if strcmp(region,'Sweden')
  slabs = slabs(1,:); % assume all the regions, had the same slabs.
end
% transmission matrix D
load Dcounties

if noNetwork
  D = sparse(zeros(size(D)));
end

test_ = test; % for plotting keep test_

% load filter data
if contains(register,'URDME') % synthetic data
  filename = mfilename('fullpath');
  Data_raw = load([filename(1:end-24) ...
                       'URDME/URDMEoutput/URDME_all']);
  Data = struct();

  if ~any(strcmp(region,regionList))
      regionList_nordic = regions(true);
      regid = strcmp(region,regionList_nordic);
      if sum(regid) == 0
          error('cannot find region in URDME data');
      end
  else
      regid = strcmp(region,regionList);
  end

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
  test = 0;
else
  Data = loadData(register);
  Data = polishData(Data,'D','Dinc',1);
  Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

  % this is to allow for perRegion posteriors
  if test >= 4
    D = 1;
    lan = 1;
    Data.D = Data.D(:, test-3);
    Data.W = Data.W(:, test-3);
    Data.H = Data.H(:, test-3);
    Npop = Npop(test-3);
    Data.regions = Data.regions(test-3);
    test = 0;
  end
end

% define filter evaluation period
ixdata = find(Data.date == datadate(1)):find(Data.date == datadate(2));
if datadate(1) < rates.meta.date(1) % data starts earlier than posterior
  warning(['data starts before posterior, assume earlier data i slab 1'...
    ', this will probably cause errors.']);
  slabs = [repmat(slabs(1),1,find(Data.date(ixdata)==rates.meta.date(1))-1)...
    slabs];
end
if numel(slabs) < numel(ixdata)
  % (posterior is from a shorter period of time)
  warning('Extending final slabs to match data.');
  slabs = [slabs repmat(slabs(end),1,numel(ixdata)-numel(slabs))];
end
ixfilter = [ixdata ixdata(end)+(1:lag)];

% change data format
Ydata = cat(3,Data.H(ixdata,:),Data.W(ixdata,:),Data.D(ixdata,:));
Ydata = permute(Ydata,[3 2 1]);

% define a common time frame, including dates
TSPAN = 1:max(ixfilter);
DATES = datenum(2e3+floor(Data.date(1)/1e4), ...
  mod(floor(Data.date(1)/1e2),1e2) ,...
  mod(Data.date(1),1e2));
DATES = DATES:DATES+numel(TSPAN)-1;
DATES = datevec(DATES);
DATES = DATES(:,1:3)*[1e4 1e2 1]'-2e7; % finally!
% selections of the above:
tspan_filter = TSPAN(ixfilter); % filter output
tspan_data = TSPAN(ixdata);  % data used
tspan_alldata = TSPAN(1:min(ixfilter(end),numel(Data.date))); % *all* data





% *** Code to allow for CSSS data, a bit measy, should be looked into.
if useCSSS
  % correction factor in xmod = fac*x
  load Fcases_interp
  fac = sum(Eave(3:end))/sum(Iave(3:end));

  % CSSS-data
  if strcmp(region,'Uppsala')
    csss = loadData('CSSS_RU');
    Ihi = csss.Ihi(:,1);
    Ilo = csss.Ilo(:,1);
    Imid = csss.Imid(:,1);
  else
    csss = loadData('CSSS');
    Ihi = csss.Ihi(:,strcmp(Data.regions, region));
    Ilo = csss.Ilo(:,strcmp(Data.regions, region));
    Imid = csss.Imid(:,strcmp(Data.regions, region));
  end
  ixcsss = find(DATES == csss.date(1)):find(DATES == csss.date(end));
  tspan_csss = TSPAN(ixcsss);

  %   Ihi = csss.Ihi;
  %   Ilo = csss.Ilo;
  %   Imid = csss.Imid;

  Imid(csss.date > 201231) = NaN; % data post 31 Dec 2020, is "bad".

  % settle on a relative variance from given asymmetric 95% CI
  varI = mean(((Ihi-Ilo)/2./Imid).^2); % generally larger
  %varI = max(varI);
  %varI = varI(test-3);
  %   varI(isnan(varI)) = [];
  %   varI = mean(varI);
  %   %varI = mean(((Ihi-Ilo)/4./Imid))^2;

  % specify observation model
  obsrates.states = {5 6 7}; % state #5, #6, #7
  obsrates.indobs = {1}; % CSSS data for the I-compartment
  obsrates.indobspars = {100./(sum(Npop)*fac)}; % probably WRONG?!
  obsrates.nstate = 8;
  obsrates.R0 = [ones(3,1); 1e-4];
  obsrates.rdiag = [(1e-3)^2*ones(3,1); varI];

  % insert CSSS data
  Idata = permute([NaN(ixcsss(1)-ixdata(1),numel(lan)); Imid; ...
    NaN(ixdata(end)-ixcsss(end),numel(lan))],[3 2 1]);
  Idata_hi = permute([NaN(ixcsss(1)-ixdata(1),numel(lan)); Ihi; ...
    NaN(ixdata(end)-ixcsss(end),numel(lan))],[3 2 1]);
  Idata_lo = permute([NaN(ixcsss(1)-ixdata(1),numel(lan)); Ilo; ...
    NaN(ixdata(end)-ixcsss(end),numel(lan))],[3 2 1]);
  Ydata = cat(1,Ydata,Idata);
else % don't use CSSS.
  % specify observation model
  obsrates.states = {5 6 7}; % state #5, #6, #7
  obsrates.indobs = {}; % No indirect measurements
  obsrates.indobspars = {};
  obsrates.nstate = 8;
  obsrates.R0 = 1*ones(3,1);
  obsrates.rdiag = (1e-3)^2*ones(3,1);
end
% specify output model
H = getC19obs(obsrates); % i.e., "H", "W", "D" ("I" if CSSS)
% read as: "Stockholm", "Uppsala", (...), "Sweden total"
%T = sparse([1 2 3*ones(1,numel(lan))],[1 2 1:numel(lan)],ones(1,2+numel(lan)));

% Add incidence compartment
rates.IStates = [1 7];
obsrates.nstate=10;
if saveall
  T = sparse(1,1:numel(lan),1);
  % I A E  PHI H W D R (Ii, Di)
  G = sparse(...
    [0 0 0 0 1 0 0 0 0 0;   % H
     0 0 0 0 0 1 0 0 0 0;   % W
     0 0 0 0 0 0 1 0 0 0;   % D
     0 0 0 0 0 0 0 0 1 0;   % Ii
     0 0 0 0 0 0 0 0 0 1;   % Di
     0 0 0 1 0 0 0 0 0 0;  % PHI
     0 0 0 0 0 0 0 1 0 0;   % R
     1 0 0 0 0 0 0 0 0 0]);   % I


  G = kron(T,G);
else
  T = [speye(numel(lan)); sparse(1,1:numel(lan),1)];
  G = kron(T,H);
end

% uncertainty parameters
Q = struct();
Q.Q0 = speye(numel(lan)*obsrates.nstate);
Q.qdiag = 0.05^2;




% these were used durring inference (some regions are sensitive)
exception.LB = -1e2;
exception.UB = 1e7;
exception.SDFAC = 0.25;
exception.AbsMagn = 1e4;

% with historic lag or not
switch type
  case 1
    % without:
    if exist([abspath(1:end-24) 'weekly/save/runs/' posterior(1:end-4) '.mat'],'file')
        warning(['Already saved prediction for ' posterior(1:end-4)]);
    end
    [Z,covZ,~,L] = C19filt([],[],G,rates,D,obsrates,Ydata,slabs, ...
                           numel(ixfilter),[],Q, exception);
    covZ.covZ = []; % save some space
                    % remove "exceptions"
    rows = find(~isinf(sum(L,1)));
    Z = Z(:,:,unique(rows));
    covZ.stdZ = covZ.stdZ(:,:,unique(rows));
    covZ.lastPriorY = covZ.lastPriorY(:,unique(rows));
    covZ.lastPostY = covZ.lastPostY(:,unique(rows));
    covZ.lastPriorYCov = covZ.lastPriorYCov(:,:,unique(rows));
    covZ.lastPostYCov = covZ.lastPostYCov(:,:,unique(rows));




    % when loading, are we loading the posterior we think we're doing?
    meta = struct();
    meta.postHash = rates.meta.hash; % same posterior
    meta.dataHash = Data.hash;  % same data
    save([abspath(1:end-24) 'weekly/save/runs/' posterior(1:end-4)],...
         'Z','covZ','meta','Ydata','tspan_filter','tspan_data','lan',...
         'useCSSS','register','DATES','TSPAN','lag','slabs','Data');
  case 2
    % with:
    if exist([abspath(1:end-24) 'weekly/save/runs/' posterior(1:end-4) '_lag.mat'],'file')
        warning(['Already saved lag prediction for ' posterior(1:end-4)]);
    end
    disp(['qdiag: ' num2str(Q.qdiag)]);
    [Z_,covZ_,Z,covZ] = ...
        C19filt_lag(G,rates,D,obsrates,Ydata,slabs,...
                    numel(ixfilter),lag,Q);
    covZ.covZ = []; % save some space
    covZ_.covZ = [];

    % when loading, are we loading the posterior we think we're doing?
    meta = struct();
    meta.postHash = rates.meta.hash; % same posterior
    meta.dataHash = Data.hash;  % same data
    save([abspath(1:end-24) 'weekly/save/runs/' posterior(1:end-4) '_lag'],...
         'Z_','covZ_','Z','covZ','meta','Ydata','tspan_filter',...
         'tspan_data','lan','useCSSS','register','DATES','TSPAN',...
         'lag','slabs','Data');
end
