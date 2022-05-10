% HORIZONSPLITCOMPARE illustrates how the dividing of optimization
% horizon is implemented. The script shares a lot of coverage with
% DYNAMIC_BETA_ML.m.
%
% The script allows for a TEST=true version, which does not generate the
% same figure as in the paper results, but instead allows for a speedy
% run-through of the code to check for potential errors. Note, the
% TEST=false versions of the script takes a long time to run, compared
% to test.


rng(0)

% if ~exist('reg','var')
%     reg=2; % region to estimate beta for
% else
%     if max(reg) > 21
%         error(['Only supported to run 21 regions,'...
%                'National posterior is sampled from a basket of region']);
%     end
% end

% if ~exist('savetofile','var')
%     savetofile = false;
% end
% if ~exist('test','var')
%     test=false;
% end
% if ~exist('verb','var')
%     verb=false;
% end

% ntime = 365;
% startdate = 200319;
% enddate = 210531;
% prevresdate = 210531; % previous result for initial guess
% if test
%     windowlength=25;
% else
%     windowlength=150;
% end
% wbeta=300000;
% steplength=20;

% % posterior
% abspath = mfilename('fullpath');
% prefix = [abspath(1:end-45) 'inference/results/KLAM/'];

% %%
% posterior = ['perRegion/slam' num2str(enddate) '_Uppsala_monthly_1_100'];
% %[postrates,~,slab] = posteriorenger(inf,[prefix posterior]);
% postrates = posteriorenger(inf,[prefix posterior]);
% slab = postrates.meta.slabs;

% % load filter data
% datasource = 'C19';
% Data = loadData(datasource);
% Data = polishData(Data,'D','Dinc',1);
% Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

% % define data and filter periods
% ixdata = find(Data.date == startdate):find(Data.date == startdate)+ntime-1;
% slab(numel(ixdata)+1:end) = [];

% ixfilter = ixdata;

% % change data format
% Ydata = cat(3,Data.H(ixdata,reg),Data.W(ixdata,reg),...
%     Data.D(ixdata,reg));
% Ydata = permute(Ydata,[3 2 1]);

% % specify observation model
% obsrates.states = {5 6 7}; % state #5, #6, #7
% obsrates.nstate = 8;
% obsrates.indobs={};
% obsrates.R0 = ones(3,1);
% obsrates.rdiag = (1e-3)^2*ones(3,1);

% G = speye(obsrates.nstate);
% [Z,covZ,~,~,S] = C19filt([],[],G,postrates,0,obsrates,Ydata,slab,numel(ixfilter));
% covZ.covZ = []; % save some space
% stdZ = covZ.stdZ; % for brevity

% nslab = max(slab);

% ntime = size(Ydata,3);
% H = getC19obs(obsrates);
% betaIdx = false(obsrates.nstate-1);
% betaIdx(3,4) = true;

% beta0 = 0.1*ones(1,ntime);


% %%
% if verb
%     display = 'iter-detailed';
% else
%     display='off';
% end
% options = optimoptions(@fmincon,'MaxFunctionEvaluations',2700000,...
%                        'MaxIterations',6000,'Display',display,...
%                        'OptimalityTolerance',2e-5);

% F=cell(nslab,1);
% for m = 1:nslab
%     C19 = getC19syst(postrates,m,1);
%     KF = getC19filt(C19);
%     F{m} = KF.F(1:7,1:7);
% end

% return;
% x0 = initC19filt(F{1},H(:,1:7),squeeze(Ydata(:,1,1)),KF.CStates(1:7));
% tic
% [x0_post,betaOpt,J,misfit,beta_mpc,ix_mpc] = mpcloop_beta_S(squeeze(Ydata),F,H(:,1:7),...
%     betaIdx,wbeta,slab,squeeze(S(:,:,:,1)),steplength,windowlength,x0,...
%     beta0,options);
% t1=toc;

if ~exist('savetofile','var')
    savetofile=false;
end

if ~exist('test','var')
    test=false;
end

if ~exist('verb','var')
    verb=1;
end

if ~exist('reg','var')
    reg=2;
end

horizon=true;

savetofile_ = savetofile;
dynamic_beta_ML;
savetofile=savetofile_;


if test
    ub = 1.05*beta0;
    lb = 0.95*beta0;
else
    ub = 5*beta0;
    lb = 0*beta0;
end
%%
tic
betaFunc = @(u)betaCost_S(u,squeeze(Ydata),F,H(:,1:7),betaIdx,wbeta,slab,squeeze(S(:,:,:,1)),x0);
[betafull,~,exitflag]=fmincon(betaFunc,beta0,[],[],[],[],...
            lb,ub,[],options);
t2=toc;
%%
if verb
    disp(['With divided horizon: ' num2str(t1)]);
    disp(['Batch optimization: ' num2str(t2)]);
end

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
    if verb
        disp(['saved figure: ' savepath2]);
    end
else
    if verb
        disp(['didn''t save figure: ' savepath2]);
    end
end
