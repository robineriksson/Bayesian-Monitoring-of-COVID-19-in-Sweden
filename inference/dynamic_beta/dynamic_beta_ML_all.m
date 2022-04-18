%DYNAMIC_BETA_ML_all Dynamic beta from posterior parameters for all regions
%   Similar to DYNAMIC_BETA_ML, but loops over all regions and collects
%   results in one file. 

% H. Runvik 2021-06-02

if ~exist('savetofile','var')
    savetofile=false;
end

%% region list
regionList = regions();

%% Options

startdate = 200401;
enddate = 210531;

wbeta=300000;
windowlength=150;
steplength=20;

% posterior
prefix = postpath;

%%
nlan=numel(regionList)

% Sweden
% Population data
load Ncounties
Weights = sum(N,1);
files = cell(size(nlan));
for r = 1:nlan
    files{r} = [prefix 'KLAM/perRegion/slam'  num2str(enddate) '_' regionList{r} '_monthly_1_100.mat'];
end

% postrates=rates;
postrates = posteriorenger(inf,files,Weights,true);
slab = postrates.meta.slabs;

% load filter data
datasource = 'C19';
Data = loadData(datasource);
Data = polishData(Data,'D','Dinc',1);
Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

% define data and filter periods
ixdata = find(Data.date == startdate):find(Data.date == enddate);
if numel(slab) < numel(ixdata)
    % (posterior is from a shorter period of time)
    warning('Extending final slab to match data.');
    slab = [slab repmat(slab(end),1,numel(ixdata)-numel(slab))];
elseif numel(slab) > numel(ixdata)
    slab(numel(ixdata)+1:end) = [];
end
ixfilter = ixdata;

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


% specify observation model
obsrates.states = {5 6 7}; % state #5, #6, #7
obsrates.nstate = 8;
obsrates.indobs={};
obsrates.R0 = ones(3,1);
obsrates.rdiag = (1e-3)^2*ones(3,1);

G = speye(obsrates.nstate);
nslab = max(slab);

H = getC19obs(obsrates);
betaIdx = false(obsrates.nstate-1);
betaIdx(3,4) = true;

options = optimoptions(@fmincon,'MaxFunctionEvaluations',1700000,...
    'MaxIterations',6000,'Display','off',...
    'OptimalityTolerance',2e-5);

slambeta=postrates.beta;

R_all = zeros(nlan,numel(ixdata));
x0_all = zeros(nlan,obsrates.nstate-1);
xSim_all = zeros(nlan,obsrates.nstate-1,numel(ixdata));

for k=1:nlan
    regionList{k}
    regiondata = cat(3,Data.H(ixdata,k),Data.W(ixdata,k),...
        Data.D(ixdata,k));
    regiondata = permute(regiondata,[3 2 1]);
    ntime = size(regiondata,3);
    postrates.beta=slambeta(:,:,k);
    [Z,covZ,~,L,S] = C19filt([],[],G,postrates,0,obsrates,regiondata,slab,numel(ixfilter));
    covZ.covZ = []; % save some space
    stdZ = covZ.stdZ; % for brevity

    beta0 = 0.1*ones(1,ntime);

    F=cell(nslab,1);
    for m = 1:nslab
        C19 = getC19syst(postrates,m,1);
        KF = getC19filt(C19);
        F{m} = KF.F(1:7,1:7);
    end
    x0 = initC19filt(F{1},H(:,1:7),squeeze(regiondata(:,1,1)),KF.CStates(1:7));

    tic
    [x0_post,betaOpt,J,misfit,beta_mpc,ix_mpc] = mpcloop_beta_S(squeeze(regiondata),F,H(:,1:7),...
        betaIdx,wbeta,slab,squeeze(S(:,:,:,1)),steplength,windowlength,x0,...
        beta0,options);
    toc

    rates_daily = dailyrates(postrates,slab);
	rates_daily.beta=betaOpt';
	R_post = getC19R0(rates_daily);
    R_all(k,:) = R_post;
    x0_all(k,:) = x0_post;

    betaStruct.vals=betaOpt';
    betaStruct.idx = betaIdx;
    xSim_all(k,:,:) = LPVsim_slabs(F,x0_post,betaStruct,slab);
end

dates = Data.date(ixdata);

postrates.dataReg = Data.reg;
postrates.dataRev = Data.rev;
postrates.dataHash = Data.hash;

postrates.beta=slambeta;

%%
savename=['dynOptPosterior' num2str(enddate) '_all_.mat'];
savename = [savename '.mat']
savename = [prefix 'dynOpt/' savename];
if savetofile
    save(savename,'R_all','dates','x0_all','postrates','xSim_all')
    disp(['saved: ' savename]);
else
    disp(['didn''t save: ' savename]);
end
