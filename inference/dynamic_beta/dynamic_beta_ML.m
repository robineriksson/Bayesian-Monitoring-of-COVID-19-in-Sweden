%DYNAMIC_BETA_ML Dynamic beta from posterior parameters.
%   Calculates a quasi-ML estimate of beta and R using MPCLOOP_BETA_S.
%   Covariance matrix S is calculated using C19FILT for a given posterior.
%   Results are optionally plotted and and saved to file.

% H. Runvik 2021-04-14

% change (sew)
if ~exist('waste', 'var')
  waste = false;
end

if ~exist('savetofile','var')
    savetofile=false;
end

if ~exist('reg','var')
    reg=2
else
    if max(reg) > 21
        error(['Only supported to run 21 regions,'...
               'National posterior is sampled from a basket of region']);
    end
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

if ~exist('windowlength','var')
    windowlength=150;
end

if ~exist('horizon','var')
    horizon=false;
end

if horizon
    ntime_ = 365;
    maxfunceval = 2700000;
else
    ntime_ = 1;
    maxfunceval = 1700000;
end

if horizon & numel(reg) > 1
    error('err:horizon','only run horizon for 1 region');
end

if ~exist('verb','var')
    verb=false;
end

if ~exist('enddate','var')
    enddate = 210531;
end

if ~exist('prevresdate','var')
    prevresdate = 210531; % previous result for initial guess
end

if ~exist('interp','var')
    interp=2;
end

if interp == 1
    ending = '_100';
elseif interp == 2
    ending ='_100_update';
else
    error('Only supporting interp 1 or 2')
end

%% region list
regionList = regions(false);

for region = reg % region to estimate beta for
    ntime = ntime_;
    if verb
        disp(['running: ' regionList{region}]);
    end

    startdate = 200401;
    if horizon
        enddate_ = startdate;
    else
        enddate_ = enddate;
    end

    wbeta=300000;
    steplength=20;

    % posterior
    prefix = postpath();

    %%

    if urdme > 0
        posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' ...
                     regionList{region} '_monthly_1_URDME' num2str(urdme) ending];
    elseif contains(register,'URDMEsew')
        if waste
          posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' ...
                       regionList{region} '_monthly_1_waste_URDMEsew_100'];
        else
          posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' ...
                       regionList{region} '_monthly_1_URDMEsew_100'];
        end
    elseif waste % change (sew)
        posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' regionList{region} '_monthly_1_waste_100'];
    else
        posterior = ['KLAM/perRegion/slam' num2str(enddate) '_' regionList{region} '_monthly_1' ending];
    end
    postrates = posteriorenger(inf,[prefix posterior]);
    slab = postrates.meta.slabs;

    % load filter data
    if contains(register,'URDMEsew') % synthetic data for
                                      % wastewater run
        filename = mfilename('fullpath');
        load([filename(1:end-38) 'sewage/urdme/data/sim_data_selected']);

        Data = struct();

        Data.regions = umod.reg;
        Data.date = umod.DATES;

        Data.H(:,region) = umod.U(5,:)';
        Data.W(:,region) = umod.U(6,:)';
        Data.D(:,region) = umod.U(7,:)';

        Data.hash = fsetop('check',Data.D(:));

        Data.rev = umod.DATES(end);
        Data.reg = register;
    elseif contains(register,'URDME') % synthetic data
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
    ixdata = find(Data.date == startdate):find(Data.date == enddate_)+ntime-1;
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

    % % define a common time frame, including dates
    % TSPAN = 1:max(ixfilter);
    % DATES = datenum(2e3+floor(Data.date(1)/1e4), ...
    %                 mod(floor(Data.date(1)/1e2),1e2) ,...
    %                 mod(Data.date(1),1e2));
    % DATES = DATES:DATES+numel(TSPAN)-1;
    % DATES = datevec(DATES);
    % DATES = DATES(:,1:3)*[1e4 1e2 1]'-2e7; % finally!
    %                                        % selections of the above:
    % tspan_filter = TSPAN(ixfilter); % filter output
    % tspan_data = TSPAN(ixdata);  % data used

    if waste
      if contains(register,'URDMEsew')
        load([filename(1:end-38) ...
              'sewage/urdme/data/simulated_sew']);
        phi_sew = WW_sim.phi_sew;
        sigma_1_period1 = WW_sim.sigma_1_period1^2;
        sigma_2_period1 = WW_sim.sigma_2_period1^2;
        sigma_1_period2 = WW_sim.sigma_1_period2^2;
        sigma_2_period2 = WW_sim.sigma_2_period2^2;

        phi_sew_data = phi_sew(ixdata);
        phi_sew_data = permute(phi_sew_data',[3 2 1]);
      else % real wastewater data:
        WW = loadData('WW');
        ixRegion = 2; % only Uppsala for now
        phi_sew = WW.phi(:,ixRegion);

        % align phi_sew with rest of data, move to function
        ixWW_start_index = find(Data.date == WW.date(1)); %
                                                          % s.t. Data.date(ixWW_start_index)
                                                          % is the
                                                          % date of
                                                          % first WW
                                                          % measurement.

        %index of last phi_sew that lies within data period
        phi_sew_stop = find(WW.date <= enddate, 1, 'last');
        %s.t. WW.date(phi_sew_stop) gives date of last WW measurement
        %within period

        ixWW_end_index = find(Data.date == WW.date(phi_sew_stop));
        %s.t. Data.date(ixWW_end_index) gives date of last WW
        %measurement within period
        ixWW_start = find(ixdata == ixWW_start_index);
        %s.t. ixdata(ixWW_start) gives starting point for WW on
        %ixdata scale
        ixWW_end = find(ixdata <= ixWW_end_index, 1, 'last');

        phi_sew_final = NaN(length(ixdata),1);
        phi_sew_final(ixWW_start:7:ixWW_end) = phi_sew(1:phi_sew_stop);

        % Append phi_sew data
        phi_sew_data = permute(phi_sew_final,[3 2 1]);
      end
        Ydata = cat(1,Ydata, phi_sew_data);

        obsrates.states = {5 6 7};
        obsrates.indobs = {4};
        rr = priorenger(inf,true,1);
        obsrates.indobspars = {rr.k_sew};
        obsrates.nstate = 8;
        load_phi_sew_noise_constants; % load sigma_1, sigma_2
                                      % Append sigma_1 and sigma_2, see sew_plots.m
        obsrates.rdiag = [(1e-3)^2*ones(3,1); sigma_2_period1^2];
        obsrates.R0 = [1*ones(3,1);sigma_1_period1^2];

        %after experimental method change:
        obsrates.rdiag_p2 = [(1e-3)^2*ones(3,1); sigma_2_period2^2];
        obsrates.R0_p2 = [1*ones(3,1);sigma_1_period2^2];
    else % don't use sewage data.

      % specify observation model
      obsrates.states = {5 6 7}; % state #5, #6, #7
      obsrates.nstate = 8;
      obsrates.indobs={};
      obsrates.R0 = ones(3,1);
      obsrates.rdiag = (1e-3)^2*ones(3,1);
    end

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

    if horizon
        beta0 = 0.1*ones(1,ntime);
    else
        try
            dyname = [prefix 'dynOpt/dynOptPosterior' num2str(prevresdate) '_' regionList{region}];
            if interp == 2
                dyname = [dyname '_update'];
            end
            oldres = load(dyname);
            oldrates_daily = dailyrates(oldres.postrates,oldres.postrates.meta.slabs');
            if isfield(oldrates_daily.meta,'interp')
                beta0 = getC19beta(oldrates_daily,oldres.R_post,oldrates_daily.meta.interp);
            else
                beta0 = getC19beta(oldrates_daily,oldres.R_post,interp);
            end
            beta0 = [beta0' beta0(end)*ones(1,ntime-numel(beta0))];
            if verb
                disp('found old guess')
            end
        catch
            beta0 = 0.1*ones(1,ntime);
            if verb
                disp('from scratch')
            end
        end
    end

    %%
    if verb
        display = 'iter-detailed';
    else
        display = 'off';
    end
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',maxfunceval,...
                           'MaxIterations',6000,'Display',display,...
                           'OptimalityTolerance',2e-5);


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
    t1 = toc;
    if verb
        t1
    end
    %%
    if horizon
        return;
    else
        rates_daily = dailyrates(postrates,slab);
        rates_daily.beta=betaOpt';
        if isfield(rates_daily.meta,'interp')
            R_post = getC19R0(rates_daily,[],rates_daily.meta.interp);
        else
            R_post = getC19R0(rates_daily,[],interp);
        end

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
            if verb
                disp('LPsim_slabs is not really working');
            end
        end

        %%
        regionList_notnordic = regions(false);
        savename=['dynOptPosterior' num2str(enddate) '_' regionList{region}];
        if contains(register,'URDMEsew')
          savename = [savename '_' register];
        elseif contains(register,'URDME')
            savename = [savename '_' register num2str(urdme)];
        end
        if waste
          savename = [savename '_waste'];
        end


        if interp == 2
            savename = [prefix 'dynOpt/' savename '_update'];
        else
            savename = [prefix 'dynOpt/' savename];
        end
        savename = [savename '.mat'];

        if savetofile
            save(savename,'R_post','dates','x0_post','postrates','xSim')
            if verb
                disp(['saved: ' savename]);
            end
        else
            if verb
                disp(['didn''t save: ' savename]);
            end
        end

    end
end
