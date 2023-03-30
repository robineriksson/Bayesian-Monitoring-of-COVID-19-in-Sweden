% PRIORPRED generates and illustrates a proxy for the prior predictive
% distribution of the Kalman filter model. This visualized "test" is
% supposed to be interpret as "could the data possibly be captured by
% the prior"?
%
% if REGEN = true, new samples are generated and save to file. If
% false, we assume that the file is available and we load it from
% file. After generating or loading, the figure is generated. See
% LAGPLOT for the posterior version of this figure.

% R. Marin 2022-04-18
rng(0)

if ~exist('Nprior','var')
    Nprior = 1e3;
end
if ~exist('savetofile','var')
    savetofile = false;
end
if ~exist('regen','var')
    regen=false;
end
if ~exist('reg','var')
    reg = 1;
else
    if max(reg) > 21
        error(['Only supported to run 21 regions,'...
               'National posterior is sampled from a basket of region']);
    end
end
if ~exist('verb','var')
    verb=false;
end

smoothD = true;
plotfun = @plot;

regionList = regions(false);
region = regionList{reg};

abspath = mfilename('fullpath');
if ~regen
    if verb
        disp(['loading priorpred_' region '.mat'])
    end
    try
        load([abspath(1:end-9) 'priorpred_' region '.mat'])
    catch
        error('file seems to be missing, run again with regen=true');
    end
else
    if verb
        disp(['generating priorpred_' region '.mat']);
    end

    lag = 7;

    useCSSS=false;    % use CSSS?
    posteriordate = '210531'; % posterior file date
    ending         = '1_100'; % at what date does the slabs start.

    % dates that "date" the data used
    datadate     = [200401 210531]; % change end date here
    register = 'C19';
    saveall      = true;
    noNetwork    = true;

    % folder for posteriors
    prefix = [postpath 'KLAM/'];

    % Population data
    load Ncounties
    Npop = sum(N,1);


    % convert region into "test" for data extraction, later.
    test = find(strcmp(region,regionList)) + 3;

    posterior = 'perRegion/';
    posterior = [posterior 'slam' posteriordate '_' region '_monthly'];
    posterior = [posterior '_' ending '.mat'];

    % sample rates
    rng(0);
    rates = posteriorenger([],[prefix posterior]);

    try
        slabs = rates.meta.slab';
    catch
        slabs = rates.meta.slabs;
    end

    rates = priorenger(Nprior,true,numel(rates.meta.slabstop)-1);
    rates.meta.date = datadate;

    % transmission matrix D
    load Dcounties


    test_ = test; % for plotting keep test_

    % load filter data

    Data = loadData(register);
    Data = polishData(Data,'D','Dinc',1);
    Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

    % this is to allow for perRegion posteriors
    if test >= 4
        D = 1;
        lan = 1;
        Data.D = Data.D(:,test-3);
        Data.W = Data.W(:,test-3);
        Data.H = Data.H(:,test-3);
        Npop = Npop(test-3);
        Data.regions = Data.regions(test-3);
        test = 0;
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


    % specify observation model
    obsrates.states = {5 6 7}; % state #5, #6, #7
    obsrates.indobs = {}; % No indirect measurements
    obsrates.indobspars = {};
    obsrates.nstate = 8;
    obsrates.R0 = 1*ones(3,1);
    obsrates.rdiag = (1e-3)^2*ones(3,1);

    % specify output model
    H = getC19obs(obsrates); % i.e., "H", "W", "D" ("I" if CSSS)
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

    % with historic lag
    [Z_,covZ_,Z,covZ] = ...
        C19filt_lag(G,rates,D,obsrates,Ydata,slabs,...
                    numel(ixfilter),lag,Q);
    covZ.covZ = []; % save some space
    covZ_.covZ = [];

    % when loading, are we loading the posterior we think we're doing?
    meta = struct();
    meta.postHash = rates.meta.hash; % same posterior
    meta.dataHash = Data.hash;  % same data

    save([abspath(1:end-9),'priorpred_' region '.mat'],...
         'Z_','covZ_','Z','covZ','meta','Ydata','tspan_filter',...
         'tspan_data','lan','useCSSS','register','DATES','TSPAN',...
         'lag','slabs','Data');
    if verb
        disp('saved priorpred.mat')
    end
end



%%
figure(1), clf;
stdZ = covZ.stdZ; % for brevity

% ***************************
% *** Postprocessing code ***
% ***************************
% 1:H 2:W 3:D 4:Ii 5:Di 6:PHI 7:R 8:I
Hcomp = 1; Wcomp = 2; Dcomp = 3; Dicomp = 5;

Icomp= [];

% Full simulation (Sweden) or single region (no transport).
Ydata_swe = squeeze(sum(Ydata,2));

% summed counties
end_ = size(Z,1);

ind = [Hcomp, Wcomp, Dicomp];

Y_swe = Z(ind,:,:);
stdY_swe = stdZ(ind,:,:);

% scaling for plots
Ydata_swe_plot = Ydata_swe;
smoothing = @(x) sgolayfilt(x,0,7);
Ydata_swe_plot(Dcomp,:) = Data.Dinc(tspan_data,find(strcmp(regionList,region)));
if smoothD
    Ydata_swe_plot(Dcomp,:) = smoothing(Ydata_swe_plot(Dcomp,:));
end

Y_swe_plot = Y_swe;


stdY_swe_plot = stdY_swe;

filtval = 1e-1;
Y_sweF = max(Y_swe,filtval); % filtered values
Y_sweF_plot = max(Y_swe_plot,filtval); % filtered values



% number of samples to draw (illustration only)
Nsamples_ = 0;


% define colors for the compartments
color = [[0 0 1]; % H
         [1 0 0]; % W
         [160 5 30]/256; % D
         [0 1 0]]; % I

alpha = 0.1;
if Nsamples_ > 0
    yyaxis left
    plotfun(tspan_filter,squeeze(Y_swe_plot(Hcomp,:,1:Nsamples_)),'-','Color',[color(1,:) alpha], ...
            'HandleVisibility','off')
    hold on;
    yyaxis right
    plotfun(tspan_filter,squeeze(Y_swe_plot(Wcomp,:,1:Nsamples_)),'-','Color',[color(2,:) alpha], ...
            'HandleVisibility','off')
    plotfun(tspan_filter,squeeze(Y_swe_plot(Dcomp,:,1:Nsamples_)),'-','Color',[color(3,:) alpha], ...
            'HandleVisibility','off')
end

% sample predictive mean

comp_vec = [Hcomp Wcomp Dcomp Icomp];

for i = comp_vec
    % note: lognormal interpretation
    if any(i == [Hcomp Icomp])
        yyaxis left, hold on;
    else
        yyaxis right, hold on;
    end
    if i == Dcomp
        plotfun(tspan_filter,mean(squeeze(Y_sweF_plot(i,:,:)),2),'-','Color',color(i,:));
    else

        y = mean(log(squeeze(Y_sweF_plot(i,:,:))),2);
        plotfun(tspan_filter,exp(y),'-','Color',color(i,:));
    end
end


% lognormal assumption: _both_ Kalman uncertainty _and_ prior
% uncertainty:

% first transform (Y,stdY) = (mean,std) interpreted as a log-normal
% variable into a normal variable (muhat,sigmahat) using the
% formulas (help lognrnd)
%   mean = exp(muhat + sigmahat^2/2)
%   std = mean * sqrt(exp(sigmahat^2) - 1)
sigmahat = sqrt(mean(log1p((stdY_swe./Y_sweF).^2),3));
sigmahat_plot = sqrt(mean(log1p((stdY_swe_plot./Y_sweF_plot).^2),3));
% (muhat,sigmahat) is then the normal transform of the Kalman filtered
% result

% then compute the prior uncertainty and add both these sources
muhat = mean(log(Y_sweF),3);
muhat_plot = mean(log(Y_sweF_plot),3);
%sigmahat = sqrt(0*sigmahat.^2+var(log(Y_sweF),0,3));
sigmahat_plot_ = sigmahat_plot;
sigmahat = sqrt(sigmahat.^2+var(log(Y_sweF),0,3));
sigmahat_plot = sqrt(sigmahat_plot.^2+var(log(Y_sweF_plot),0,3));



kappa = 1;
Ylo = exp(muhat-kappa*sigmahat);
Yhi = exp(muhat+kappa*sigmahat);

Ylo_plot = exp(muhat_plot-kappa*sigmahat_plot);
Yhi_plot = exp(muhat_plot+kappa*sigmahat_plot);

% minor filtering to avoid issues in log-scale and patch
Ylo = real(Ylo); Ylo(Ylo < 0) = 0;
Yhi = real(Yhi); Yhi(Yhi < 0) = 0;

Ylo_plot = real(Ylo_plot); Ylo_plot(Ylo_plot < 0) = 0;
Yhi_plot = real(Yhi_plot); Yhi_plot(Yhi_plot < 0) = 0;

tt = [tspan_filter fliplr(tspan_filter)];


for i = comp_vec
    % note: lognormal interpretation
    if any(i == [Hcomp Icomp])
        yyaxis left
    else
        yyaxis right
    end
    %patch(tt,max(1,[Ylo_plot(i,:) fliplr(Yhi_plot(i,:))]), ...
    patch(tt,max(0,[Ylo_plot(i,:) fliplr(Yhi_plot(i,:))]), ...
          [0.9 0.9 0.9],'FaceAlpha',0.15, ...
          'LineStyle','none','FaceColor',color(i,:),...
          'HandleVisibility','off');
end



% data used in filter
yyaxis left
plotfun(tspan_data,Ydata_swe_plot(1,:),'.','Color',color(1,:),'HandleVisibility','off')
yyaxis right
plotfun(tspan_data,Ydata_swe_plot(2,:),'.','Color',color(2,:),'HandleVisibility','off')
yyaxis right
plotfun(tspan_data,Ydata_swe_plot(3,:),'.','Color',color(3,:),'HandleVisibility','off')


% slab limiters
slabstops = find(diff(slabs))+1;
if strcmp(region,'Sweden')
    slabstop = slabstops(:,1);
end

for k = slabstops
    xline(tspan_data(k),':k');
end
hold off;

xlim(tspan_data([1 end]))
xtk = fliplr(tspan_data(end:-28:1));

if reg == 2
    yyaxis right % no negative axis
    ylim([0 40])

    yyaxis left
    ylim([0,200])
elseif reg == 1
        yyaxis left
        ylim([0,1200])
        yyaxis right
        ylim([0, 220])
end




set(gca,'xtick',xtk);


% wanted xticks
% for abbreviations, see https://www.bydewey.com/monthdayabb.html
strmonths = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'June'...
             'July' 'Aug' 'Sept' 'Oct' 'Nov' 'Dec'};

% 15th each month:
ixdates = find(mod(DATES(tspan_filter),100) == 1);
xxdates = tspan_filter(ixdates);
dates=DATES(xxdates);
xticks(xxdates);
% what month:
months = mod(floor(dates/1e2),1e2);
slabtitle = strmonths(months);
slabtitle(2:2:end) = {''}; % only every 2nd

% xticks & -labels
xticklabels(slabtitle);
xtickangle(45);


% 2021 delimiter
xline(find(DATES(tspan_filter) == 210101)+tspan_filter(1)-1,'--k','2021', ...
      'LabelVerticalAlignment','top','LabelOrientation','horizontal', ...
      'interpreter','latex','HandleVisibility','off');



leg = {'Hospital [H]','Intensive [W] $\,$(right)','Dead/day [D] (right)'};

numcols = 3;


legend(leg, ...
       'Location','northwest','interpreter','latex','FontSize',12,...
       'Orientation','horizontal',...
       'NumColumns', 1);




hold off;



savepath = abspath(1:end-17);
figname = [savepath 'priorpred_' region '.pdf'];
if savetofile % true if we should compute the error table



    % finalize the plot
    h = gcf;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';

    % polish target output size
    set(h,'PaperPositionMode','auto');
    set(h,'Position',[100 100 500 350]);

    %print(h, ['~/Desktop/lag_' region '.eps'], '-depsc')
    print(h, figname, '-dpdf')
    if verb
        disp(['saved figure:' figname])
    end
else
    if verb
        disp(['didn''t save figure:' figname])
    end
end