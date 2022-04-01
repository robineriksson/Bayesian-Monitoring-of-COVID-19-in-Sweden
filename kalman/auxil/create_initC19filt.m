function [Y0, Y0Cov, meta] = create_initC19filt(Nreplicas,file,date)
%CREATE_INITC19FILT generates a savepoint (state posterior including hidden states)
%  that allows one to start from a later state, i.e., instead of restarting
%  from data0 this creates a new zeropoint.
%
%   [Y0,Y0COV] = CREATE_INITC19FILT(Nreplicas,FILE,date) computes the mean
%   out of Nreplicas from posterior FILE. Are in the correct format to be
%   used in C19FILT(Y0,Y0Cov,...).
%
%   NREPLICAS number of rates to sample. Can be left empty {[]} and will
%   then call for all rates in the posterior from FILE, see POSTERIORENGER.
%
%   FILE is the name of the posterior .mat file. The naming should
%   preferably follow the pattern 'abcYYMMDD_XX.mat', where XX
%   corresponds to [low/mid/high, low/mid/high] and is represented by
%   [1/2/3 1/2/3]. The file should preferably be stored in
%   /inference/results/.
%
%   DATE the end date that we want to run the simulations for.
%
%   [...,meta] = CREATE_INITC19FILT(...) returns a meta file holding
%   * date   (end date in number format, i.e., 200501)
%   * region (what region, e.g., Uppsala)
%
%   Example 1:
%      % filename.
%      file = ['inference/results/SLAM/slam210331_Uppsala_'...
%               'monthly_1_100.mat'];
%      % create origo point.
%      [Y0,Y0Cov] = create_initC19filt([],file,210331);
%      % start from the origo point.
%      [Z,covZ] = C19filt(Y0,Y0Cov,...); % add missing inputs here
%

% R. Eriksson 05-12-21

% from what date to start from. Could consider running from 200401 as well.
datestart = 200319;

% use CSSS?
useCSSS=contains(file,'csss');

% sample rates
rates = posteriorenger(Nreplicas,file);
slab = rates.meta.slab';

% transmission matrix D
load Dcounties

% Population data
load Ncounties
Npop0 = sum(N,1);

% load filter data
% *** find datasource from posterior?
datasource = 'RU';
Data = loadData(datasource);
Data = polishData(Data,'D','Dinc',1);
Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

% this is to allow for perRegion posteriors
D = 1;
lan = 1;
reg = find(strcmp(Data.regions,rates.meta.region));
Data.D = Data.D(:,reg);
Data.W = Data.W(:,reg);
Data.H = Data.H(:,reg);
Npop = Npop0(reg);
Data.regions = Data.regions(reg);


% define filter evaluation period
ixdata = find(Data.date == datestart):find(Data.date == date);
Data.date = Data.date(ixdata);

if datestart(1) < rates.meta.date(1) % data starts earlier than posterior
  warning('data starts before posterior, assume earlier data i slab 1');
  slab = [repmat(slab(1),1,find(Data.date==rates.meta.date(1))-1)...
    slab];
end

if numel(slab) < numel(ixdata)
  % (posterior is from a shorter period of time)
  warning('Extending final slab to match data.');
  slab = [slab repmat(slab(end),1,numel(ixdata)-numel(slab))];
end


% change data format
Ydata = cat(3,Data.H(ixdata,:),Data.W(ixdata,:),Data.D(ixdata,:));
Ydata = permute(Ydata,[3 2 1]);


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
    Ihi = csss.Ihi(:,reg);
    Ilo = csss.Ilo(:,reg);
    Imid = csss.Imid(:,reg);
  end
  ixcsss = find(Data.date(1) == csss.date(1)):find(Data.date(end) == csss.date(end));
  
  Imid(csss.date > 201231) = NaN; % data post 31 Dec 2020, is "bad".
  
  % settle on a relative variance from given asymmetric 95% CI
  varI = mean(((Ihi-Ilo)/2./Imid).^2); % generally larger
  
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
T = [speye(numel(lan)); sparse(1,1:numel(lan),1)];
G = kron(T,H);


% uncertainty parameters
Q = struct();
Q.Q0 = speye(numel(lan)*obsrates.nstate);
Q.qdiag = 0.05^2;

% these were used durring inference (some regions are sensitive)
exception.LB = -5e2;
exception.UB = 1e7;
exception.SDFAC = 0.25;
exception.AbsMagn = 1e4;

% run filter
[Z,covZ,~,L] = C19filt([],[],G,rates,D,obsrates,Ydata,slab, ...
  numel(ixdata),[],Q,exception);
if any(isinf(L(:)))
  warning('exceptions found');
end
% mean of all the runs.
dims = size(covZ.lastPostY);
Y0 = mean(covZ.lastPostY,numel(dims));
Y0Cov = mean(covZ.lastPostYCov,numel(dims)+1); % by definition 1 dim larger

% *** meta struct: date, region
meta = struct('date',{date},'region',rates.meta.region);

