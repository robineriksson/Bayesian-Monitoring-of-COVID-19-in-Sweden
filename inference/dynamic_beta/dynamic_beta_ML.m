%DYNAMIC_BETA_ML Dynamic beta from posterior parameters

% H. Runvik 2021-04-14

if ~exist('savetofile','var')
    savetofile=false;
end

if ~exist('reg','var')
    reg=2
end

if ~exist('evalplot','var')
    evalplot=false;
end

if ~exist('urdme','var')
    urdme=0;
end

if ~exist('register','var')
    register='C19';
end

%% region list
regionList = regions(false);



for region = reg % region to estimate beta for
    disp(['running: ' regionList{region}]);
    startdate = 200401;
    enddate = 210531;
    prevresdate = 210531; % previous result for initial guess

    wbeta=300000;
    windowlength=150;
    steplength=20;

    % posterior
    prefix = postpath();

    %%

    if urdme > 0
        posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' regionList{region} '_monthly_1_URDME' num2str(urdme) '_100'];
    else
        posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' regionList{region} '_monthly_1_100'];
    end
    postrates = posteriorenger(inf,[prefix posterior]);
    slab = postrates.meta.slabs;

    % load filter data
    if contains(register,'URDME') % synthetic data
        filename = mfilename('fullpath');
        Data_raw = load([filename(1:end-38) ...
                         'URDME/URDMEoutput/URDME_all']);

        Data = struct();
        Data.regions = Data_raw.D.regions;
        Data.date = Data_raw.D.date;
        % which data "page"
        try
            runid = str2double(register(end));
        catch
            runid = 1;
        end
        Data.H = squeeze(Data_raw.D.U(:,5,:,runid));
        Data.W = squeeze(Data_raw.D.U(:,6,:,runid));
        Data.D = squeeze(Data_raw.D.U(:,7,:,runid));

        Data.hash = fsetop('check',Data_raw.D.U(:));

        Data.rev = Data_raw.D.date(end);
        Data.reg = register;
    else
        Data = loadData(register);
        Data = polishData(Data,'D','Dinc',1);
        Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});
    end
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

    % change data format
    Ydata = cat(3,Data.H(ixdata,region),Data.W(ixdata,region),...
                Data.D(ixdata,region));
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

    % specify observation model
    obsrates.states = {5 6 7}; % state #5, #6, #7
    obsrates.nstate = 8;
    obsrates.indobs={};
    obsrates.R0 = ones(3,1);
    obsrates.rdiag = (1e-3)^2*ones(3,1);

    G = speye(obsrates.nstate);
    D = 0;

    [Z,covZ,~,~,S] = C19filt([],[],G,postrates,D,obsrates,Ydata(:,1,:),slab,numel(ixfilter));
    covZ.covZ = []; % save some space
    stdZ = covZ.stdZ; % for brevity

    nslab = max(slab);

    ntime = size(Ydata,3);
    H = getC19obs(obsrates);
    betaIdx = false(obsrates.nstate-1);
    betaIdx(3,4) = true;

    try
        oldres = load([prefix 'dynOpt/dynOptPosterior' num2str(prevresdate) '_'...
                       regionList{region} '.mat']);
        oldrates_daily = dailyrates(oldres.postrates,oldres.postrates.meta.slabs');
        beta0 = getC19beta(oldrates_daily,oldres.R_post);
        beta0 = [beta0' beta0(end)*ones(1,ntime-numel(beta0))];
        disp('found old guess')
    catch
        beta0 = 0.1*ones(1,ntime);
        disp('from scratch')
    end

    %%
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',1700000,...
                           'MaxIterations',6000,'Display','off',...'iter-detailed',...
                           'OptimalityTolerance',2e-5);
    %options.Algorithm = 'sqp';

    F=cell(nslab,1);
    for m = 1:nslab
        C19 = getC19syst(postrates,m,1);
        KF = getC19filt(C19);
        F{m} = KF.F(1:7,1:7);
    end

    x0 = initC19filt(F{1},H(:,1:7),squeeze(Ydata(:,1,1)),KF.CStates(1:7));
    tic
    [x0_post,betaOpt,J,misfit,beta_mpc,ix_mpc] = mpcloop_beta_S(squeeze(Ydata(:,1,:)),F,H(:,1:7),...
                                                                betaIdx,wbeta,slab,squeeze(S(:,:,:,1)),steplength,windowlength,x0,...
                                                                beta0,options);
    toc
    %%
    rates_daily = dailyrates(postrates,slab);
    rates_daily.beta=betaOpt';
    R_post = getC19R0(rates_daily);

    dates = Data.date(ixdata);

    postrates.dataReg = Data.reg;
    postrates.dataRev = Data.rev;
    postrates.dataHash = Data.hash;

    %%

    if evalplot
        figure(1), clf
        plot(R_post)
    end

    beta_mpc(end,:)=[];
    ix_mpc(end,:)=[];

    %%



    if evalplot

        beta_temp = 0.1*ones(size(betaOpt,2),size(beta_mpc,1));
        for k=1:size(beta_mpc,1)
            beta_temp(ix_mpc(k,:),k)=beta_mpc(k,:);
        end
        rates_daily.beta=beta_temp;
        R_mpc = getC19R0(rates_daily);

        Ic_diff=zeros(size(beta_mpc));
        beta_diff=zeros(size(beta_mpc));


        for k=1:size(beta_mpc,1)
            Ic_diff(k,:) = R_mpc(ix_mpc(k,:),k)-R_post(ix_mpc(k,:));
            beta_diff(k,:) = beta_mpc(k,:)-betaOpt(ix_mpc(k,:));
        end

        figure(2), clf
        plot(Ic_diff','color',[0.9 0.9 1])
        yline(0,':');
        xline(steplength,'--');
        hold on
        plot(mean(Ic_diff),'b')
        plot(mean(Ic_diff)+std(Ic_diff),'--b')
        plot(mean(Ic_diff)-std(Ic_diff),'--b')
        hold off


        figure(3), clf
        plot(beta_diff','color',[0.9 0.9 1])
        yline(0,':');
        xline(steplength,'--');
        hold on
        plot(mean(beta_diff),'b')
        plot(mean(beta_diff)+std(beta_diff),'--b')
        plot(mean(beta_diff)-std(beta_diff),'--b')
        hold off

        figure(4), clf
        plot(ix_mpc',beta_mpc')
    end
    %%
    try
        % betaStruct.vals=betaOpt';
        % betaStruct.idx = betaIdx;
        %xSim = LPVsim_slabs(F,x0_post,betaStruct,slab);
        xSim = LPVsim_slabs_par(F,x0_post,betaOpt',betaIdx,slab);
        if evalplot
            figure(5), clf
            subplot(3,1,1)
            plot(xSim(5,:))
            hold on
            plot(squeeze(Ydata(1,1,:)))
            hold off
            subplot(3,1,2)
            plot(xSim(6,:))
            hold on
            plot(squeeze(Ydata(2,1,:)))
            hold off
            subplot(3,1,3)
            plot(xSim(7,:))
            hold on
            plot(squeeze(Ydata(3,1,:)))
            hold off
        end
    catch
        disp('LPsim_slabs is not really working');
    end

    %%
    regionList_notnordic = regions(false);
    savename=['dynOptPosterior' num2str(enddate) '_' regionList{region}];
    if urdme
        savename = [savename '_URDME' num2str(urdme)];
    end
    savename = [savename '.mat']
    savename = [prefix 'dynOpt/' savename];
    if savetofile
        save(savename,'R_post','dates','x0_post','postrates','xSim')
        disp(['saved: ' savename]);
    else
        disp(['didn''t save: ' savename]);
    end

end



return;