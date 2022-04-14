rng(0)

if ~exist('reg','var')
    reg=2; % region to estimate beta for
end
if ~exist('savetofile','var')
    savetofile = false;
end

ntime = 365;
startdate = 200319;
enddate = 210531;
prevresdate = 210531; % previous result for initial guess

wbeta=300000;
windowlength=150;
steplength=20;

% posterior
abspath = mfilename('fullpath');
prefix = [abspath(1:end-45) 'inference/results/KLAM/'];

%%
posterior = ['perRegion/slam' num2str(enddate) '_Uppsala_monthly_1_100'];
%[postrates,~,slab] = posteriorenger(inf,[prefix posterior]);
postrates = posteriorenger(inf,[prefix posterior]);
slab = postrates.meta.slabs;

% load filter data
datasource = 'C19';
Data = loadData(datasource);
Data = polishData(Data,'D','Dinc',1);
Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

% define data and filter periods
ixdata = find(Data.date == startdate):find(Data.date == startdate)+ntime-1;
slab(numel(ixdata)+1:end) = [];

ixfilter = ixdata;

% change data format
Ydata = cat(3,Data.H(ixdata,reg),Data.W(ixdata,reg),...
    Data.D(ixdata,reg));
Ydata = permute(Ydata,[3 2 1]);

% specify observation model
obsrates.states = {5 6 7}; % state #5, #6, #7
obsrates.nstate = 8;
obsrates.indobs={};
obsrates.R0 = ones(3,1);
obsrates.rdiag = (1e-3)^2*ones(3,1);

G = speye(obsrates.nstate);
[Z,covZ,~,~,S] = C19filt([],[],G,postrates,0,obsrates,Ydata,slab,numel(ixfilter));
covZ.covZ = []; % save some space
stdZ = covZ.stdZ; % for brevity

nslab = max(slab);

ntime = size(Ydata,3);
H = getC19obs(obsrates);
betaIdx = false(obsrates.nstate-1);
betaIdx(3,4) = true;

beta0 = 0.1*ones(1,ntime);


%%
options = optimoptions(@fmincon,'MaxFunctionEvaluations',2700000,...
                       'MaxIterations',6000,'Display',...
                       'iter-detailed',...                       %'off',...
                       'OptimalityTolerance',2e-5);

F=cell(nslab,1);
for m = 1:nslab
    C19 = getC19syst(postrates,m,1);
    KF = getC19filt(C19);
    F{m} = KF.F(1:7,1:7);
end

x0 = initC19filt(F{1},H(:,1:7),squeeze(Ydata(:,1,1)),KF.CStates(1:7));
tic
[x0_post,betaOpt,J,misfit,beta_mpc,ix_mpc] = mpcloop_beta_S(squeeze(Ydata),F,H(:,1:7),...
    betaIdx,wbeta,slab,squeeze(S(:,:,:,1)),steplength,windowlength,x0,...
    beta0,options);
t1=toc;

ub = 5*beta0;
lb = 0*beta0;
%%
tic
betaFunc = @(u)betaCost_S(u,squeeze(Ydata),F,H(:,1:7),betaIdx,wbeta,slab,squeeze(S(:,:,:,1)),x0);
[betafull,~,exitflag]=fmincon(betaFunc,beta0,[],[],[],[],...
            lb,ub,[],options);
t2=toc;
%%
disp(['With divided horizon: ' num2str(t1)]);
disp(['Batch optimization: ' num2str(t2)]);

%%
figure(1)
plot(betafull,'b')
hold on
plot(ix_mpc',beta_mpc','r--')
hold off
xticklabels([])
grid
ylabel('$\beta_t$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')
legend({'Batch optimization','Divided horizon'},'Interpreter','Latex')

savepath = mfilename('fullpath');
savepath2 = [savepath(1:end-27) 'horizonTest'];
if savetofile
    % finalize the plot
    h2 = gcf;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    % polish target output size
    set(h2,'PaperPositionMode','auto');
    set(h2,'Position',[100 100 500 350]);
    print('-dpdf', savepath2)
    disp(['saved figure: ' savepath2]);
else
    disp(['didn''t save figure: ' savepath2]);
end
