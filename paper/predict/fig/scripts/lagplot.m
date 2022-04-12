%WEEKLY_PLOTS generates plots for the weekly report, using generated data
%  from WEEKLY_PREDICT

% Generates two plots per region
% type 1: 14 days ahead prediciton.
% type 2: 7 days ahead historic prediction.
%
% and the regions of interests are{'Uppsala' 'Stockholm' 'Sweden'}


if ~exist('generateData','var')
    generateData = 0; % 1 if you need to generate the dataset
end
if generateData
    laggen;
end

smoothD = 1; % 1 if SG filter of Dinc

% this controls the update of .tex- and .pdf-files (so not
% .mat-files):
if ~exist('savetofile','var')
    savetofile = 0;
end

% ****************************************************
% *** Change values here for different loaded posterior predictions
% ****************************************************


% posterior data
posteriordate = '210531';
ending        = '1_100'; % at what date does the slabs start.



Iplotting   = false; % include I compartment in plots
                     %phiplotting = false; % additional plot for phi compartment
                     % ****************************************************
                     % *** Don't change below, if you only want other plots.
                     % ****************************************************

% plot-function: @plot
plotfun = @plot;

if ~exist('reg','var')
    reg = [2];
end
regions_ = regions();
j = 0;
for rid = reg
    region = regions_(rid)
    j = j + 1;
    region = region{:} % to character

    figure(j), clf;

    % construct posterior file name
    % folder for posteriors
    abspath = mfilename('fullpath');
    prefix = [postpath '/SLAM/'];
    if ~strcmp(region,'Sweden')
        posterior = 'perRegion/';
        lan = {1};
    else
        posterior = '';
    end

    posterior = [posterior 'slam' posteriordate '_' region '_monthly_' ...
                 ending];
    rates = posteriorenger([],[prefix posterior]);



    try
        load([abspath(1:end-33) 'weekly/save/runs/' posterior '_lag']);
    catch % missing file. Recompute?
           code =  ['reg=' num2str(rid) ', ', ...
                 'type=2, '...
                 'laggen'];
        error('MATLAB:ambiguousSyntax',['does not find prediction file. Did you run weekly_predict?'...
                                        '\nTry: \n' code]);
    end

    if any(meta.postHash ~= rates.meta.hash)
        warning('the saved data do not match the data or posterior');
    end


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
    Ydata_swe_plot(Dcomp,:) = Data.Dinc(tspan_data,find(strcmp(regions,region)));
    if smoothD
        Ydata_swe_plot(Dcomp,:) = smoothing(Ydata_swe_plot(Dcomp,:));
    end

    Y_swe_plot = Y_swe;
    stdY_swe_plot = stdY_swe;

    filtval = 0;
    Y_sweF = max(Y_swe,filtval); % filtered values
    Y_sweF_plot = max(Y_swe_plot,filtval); % filtered values



    % number of samples to draw (illustration only)
    Nsamples_ = 10;

    % define colors for the compartments
    color = [[0 0 1]; % H
             [1 0 0]; % W
             [160 5 30]/256; % D
             [0 1 0]]; % I

    yyaxis left
    plotfun(tspan_filter,squeeze(Y_swe_plot(Hcomp,:,1:Nsamples_)),'-','Color',[color(1,:) 0.1], ...
            'HandleVisibility','off')
    hold on,
    yyaxis right
    plotfun(tspan_filter,squeeze(Y_swe_plot(Wcomp,:,1:Nsamples_)),'-','Color',[color(2,:) 0.1], ...
            'HandleVisibility','off')
    if type == 2
        plotfun(tspan_filter,squeeze(Y_swe_plot(Dcomp,:,1:Nsamples_)),'-','Color',[color(3,:) 0.1], ...
                'HandleVisibility','off')
        if ~isempty(Icomp)
            plotfun(tspan_filter,squeeze(Y_swe_plot(Icomp,:,1:Nsamples_)),'-','Color',[color(4,:) 0.1], ...
                    'HandleVisibility','off')
        end
    end

    % sample predictive mean
    if type == 2
        comp_vec = [Hcomp Wcomp Dcomp Icomp];
    else
        comp_vec = [Hcomp Wcomp];
    end
    for i = comp_vec
        % note: lognormal interpretation
        if any(i == [Hcomp Icomp])
            yyaxis left
        else
            yyaxis right
        end
        plotfun(tspan_filter,exp(mean(log(squeeze(Y_sweF_plot(i,:,:))),2)),'-','Color',color(i,:));
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
    Ylo = exp(muhat-sigmahat);
    Yhi = exp(muhat+sigmahat);

    Ylo_plot = exp(muhat_plot-sigmahat_plot);
    Yhi_plot = exp(muhat_plot+sigmahat_plot);


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
    if type == 2
        yyaxis right
        plotfun(tspan_data,Ydata_swe_plot(3,:),'.','Color',color(3,:),'HandleVisibility','off')
        if size(Ydata_swe_plot,1) == 4
            yyaxis left
            plotfun(tspan_data,Ydata_swe_plot(4,:),'.','Color',color(4,:),'HandleVisibility','off')
        end
    end

    % slab limiters
    slabstops = find(diff(slabs))+1;
    if strcmp(region,'Sweden')
        slabstop = slabstops(:,1);
    end

    for k = slabstops
        xline(tspan_data(k),':k');
    end
    hold off,


    xlim(tspan_data([1 end]))
    xtk = fliplr(tspan_data(end:-28:1));
    yyaxis right % no negative axis
    yvr = ylim;
    ylim([0 yvr(2)])





    % below is a trick to make the grid align both with the left and right
    % yaxis.

    yyaxis left
    ytl = get(gca, 'YTick'); % Get Controlling Left Ticks
    ytlv = round(linspace(min(ylim),max(ylim),numel(ytl)));
    ytlc = compose('%d',ytlv); % Tick Label Cell Array
    set(gca, 'YTick',ytl, 'YTickLabel',ytlc)  % Set New right Tick Labels



    yyaxis right
    ylim0 = ylim;
    ytr = get(gca, 'YTick'); % Get Right Tick Values

    tinc = round(diff(ylim)/numel(ytl));  % Tick Increment


    yyaxis right
    tinc0 = round(diff(ylim)/numel(ytl));  % Tick Increment
    ytrv = min(ylim):tinc0:(min(ylim) + tinc0*(numel(ytl)-1));



    ylim([min(ytrv) max(ytrv)])% Set New ‘ylim’
    ytrc = compose('%d',ytrv); % Tick Label Cell Array
    set(gca, 'YTick',ytrv, 'YTickLabel',ytrc)  % Set New Right Tick Labels
    grid on

    if min(ylim) + tinc0*(numel(ytl)-1) < max(ylim0)
        % top is missing, add one tick on left and right.


        yyaxis left
        ytl = get(gca, 'YTick'); % Get Controlling Left Ticks
        tinc = diff(ytl(1:2));
        ytlv = [ytl ytl(end) + tinc];
        ytlc = compose('%d',ytlv); % Tick Label Cell Array
        set(gca, 'YLim', [min(ytlv) max(ytlv)])
        set(gca, 'YTick',ytlv, 'YTickLabel',ytlc)  % Set New left Tick Labels


        yyaxis right
        ytr = get(gca, 'YTick'); % Get Controlling Left Ticks
        tinc = diff(ytr(1:2));
        ytrv = [ytr ytr(end) + tinc];
        ytrc = compose('%d',ytrv); % Tick Label Cell Array
        set(gca, 'YLim', [min(ytrv) max(ytrv)])
        set(gca, 'YTick',ytrv, 'YTickLabel',ytrc)  % Set New Right Tick Labels   end
    end
    % end of trick.


    %% manually changing the limits!
    if rid == 2
        yyaxis left
        ylim([0,160])
        yyaxis right
        ylim([0, 32])
    else
        yyaxis left
        ylim([0,1200])
        yyaxis right
        ylim([0, 220])
    end

    %%

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




    % ***************************************
    % *** Create tables and save to files ***
    % ***************************************
    savepath = mfilename('fullpath');
    savepath = [savepath(1:end-15) 'lag_' region];
    if savetofile % true if we should compute the error table



        % finalize the plot
        h = gcf;
        ax = gca;
        ax.TickLabelInterpreter = 'latex';

        % polish target output size
        set(h,'PaperPositionMode','auto');
        set(h,'Position',[100 100 500 350]);

        print('-dpdf', savepath)
        disp(['saved figure: ' savepath]);


    end
end
