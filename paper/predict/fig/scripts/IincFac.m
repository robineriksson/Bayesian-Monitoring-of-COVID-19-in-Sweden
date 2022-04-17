%RECPLOT generates the (R)ecovered estimate plot and TABLE
%
% How to use:
% 1) run weekly_predict for:
%    a) Uppsala, Stockholm, Sweden
%    b) type 1 (lag 14 ahead)
%    c) specific posterior and dates
% 2) Run this script, and the final table, table_R is in latex form,
%    or have a look at datastore.table_ci for the credible interval in cell
%    string format.
%

% R. Eriksson 2021-05-12.

% ****************************************************
% *** Change values here for different loaded posterior predictions
% ****************************************************

% posterior data
posteriordate = '210531';
ending        = '1_100'; % at what date does the slabs start.

% this controls the update of .tex- and .pdf-files (so not
% .mat-files):
if ~exist('savetofile','var')
    savetofile = false;
end


% which regions to include (see regionList below)
if ~exist('reg','var')
    reg = [1 2];
end


% ****************************************************


regionList = regions(false);

datastore = struct();


for rid = reg
    region = regionList{rid}; % to character

    % construct posterior file name
    % folder for posteriors
    abspath = mfilename('fullpath');
    prefix = [postpath 'KLAM/'];
    if ~strcmp(region,'Sweden')
        posterior = 'perRegion/';
        lan = {1};
    else
        posterior = '';
    end

    posterior = [posterior 'slam' posteriordate '_' region '_monthly_' ...
                 ending];
    rates = posteriorenger([],[prefix posterior]);


    % load weekly prediction.
    pp = postpath();
    try
        load([pp(1:end-18) 'weekly/save/runs/' posterior]);
    catch % missing file. Recompute?
        code =  ['reg=' num2str(rid) ', ', ...
                 'type=1, '...
                 'laggen'];
        error('MATLAB:ambiguousSyntax',['does not find prediction file. Did you run weekly_predict?'...
                                        '\nTry: \n' code]);
    end

    if any(meta.postHash ~= rates.meta.hash)
        error('the saved data do not match the data or posterior');
    end

    % Separating R from the rest of Z
    ixI = 4;
    ZI = Z(ixI,:,:);
    stdI = covZ.stdZ(ixI,:,:);



    %% Compute (R)ecovered.

    % Population data
    load Ncounties
    Npop = sum(N,1);
    regi = find(ismember(regionList,region));


    % Uncertainty
    stdYI = stdI;
    YI = ZI;
    sigmahat_I = sqrt(mean(log1p((stdYI./YI).^2),3));
    muhat_I = mean(log(YI),3);

    sigmahat_tot = sqrt(sigmahat_I.^2+var(log(YI),0,3));


    datastore.(region).YR = YI;
    datastore.(region).muhat = muhat_I;
    datastore.(region).sigmahat = sigmahat_tot;
    datastore.(region).DATES = DATES;
    datastore.(region).tspan_filter = tspan_filter;

    %% process for table
    datastore.(region).low95 = exp(datastore.(region).muhat - 1.96*datastore.(region).sigmahat);
    datastore.(region).low68 = exp(datastore.(region).muhat - 1.0*datastore.(region).sigmahat);
    datastore.(region).mid = exp(datastore.(region).muhat);
    datastore.(region).high68 = exp(datastore.(region).muhat + 1.0*datastore.(region).sigmahat);
    datastore.(region).high95 = exp(datastore.(region).muhat + 1.96*datastore.(region).sigmahat);



    dates = DATES(tspan_data);
    t = 1:numel(DATES);



end

%% Plot and compare incidence with (unseen, but correlated) data
% load and compare data with FHM
data_inc = loadData('FHM_Iinc'); % number of positive tests (incidence)
data_test = loadData('FHM_test'); % number of tests (total tests)
%%
% Policy changes regarding testing
% see https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/statistik-och-analyser/bekraftade-fall-i-sverige/
policy_dates = {200313 200605 200831 211101 211122 220209};


for ri = reg
    region = regionList{ri}; % to character

    DATES = datastore.(region).DATES;
    t = 1:numel(DATES);

    %   if numel(reg) > 1
    %     disp([num2str(ri) ': ' region '. ' '<key>']);
    %     pause;
    %   end
    %


    figure(1), clf,  hold on

    data_c = [1 0 0];
    filter_c = [0 0 1];
    day_smooth = 7;
    smoothing = @(x,k) sgolayfilt(x,0,k);

    % Illustrate the testing data
    if data_inc.date(end) > DATES(end)
        tspan_fhm = find(DATES == data_inc.date(1)):numel(DATES);
    else
        tspan_fhm = find(DATES == data_inc.date(1)):find(DATES == data_inc.date(end));
    end

    ydata = smoothing(data_inc.Iinc(tspan_fhm,ri),day_smooth);
    ystd = movstd(data_inc.Iinc(tspan_fhm,ri),day_smooth);

    % select only the data we need for now.
    plot(t(tspan_fhm),ydata,'-','Color',data_c, 'HandleVisibility','on'); hold on
    s=scatter(t(tspan_fhm),data_inc.Iinc(tspan_fhm,ri),'filled', 'MarkerFaceAlpha',0.2,...
              'MarkerFaceColor',data_c,'SizeData',10,'HandleVisibility','off');

    data_68 = [[ydata-ystd]',fliplr([ydata+ystd]')];
    data_95 = [[ydata-2*ystd]',fliplr([ydata+2*ystd]')];


    pred_tt =  [t(tspan_fhm), fliplr(t(tspan_fhm))];
    patch(pred_tt, data_68,...
          [1 1 1],'FaceAlpha',0.15, ...
          'LineStyle','none','FaceColor',data_c, 'HandleVisibility','off');
    patch(pred_tt, data_95,...
          [1 1 1],'FaceAlpha',0.05, ...
          'LineStyle','none','FaceColor',data_c, 'HandleVisibility','off');


    % illustrate the
    tspan_filter = datastore.(region).tspan_filter;
    plot(t(tspan_filter),datastore.(region).mid,'Color',filter_c)
    pred_68 = [datastore.(region).low68, fliplr(datastore.(region).high68)];
    pred_95 = [datastore.(region).low95, fliplr(datastore.(region).high95)];

    pred_tt =  [t(tspan_filter), fliplr(t(tspan_filter))];
    patch(pred_tt, pred_68,...
          [1 1 1],'FaceAlpha',0.15, ...
          'LineStyle','none','FaceColor',filter_c, 'HandleVisibility','off');

    patch(pred_tt,pred_95, ...
          [1 1 1],'FaceAlpha',0.05, ...
          'LineStyle','none','FaceColor',filter_c, 'HandleVisibility','off');





    % Testing policy change
    for i = 1:numel(policy_dates)
        if policy_dates{i} < DATES(end)
            xline(find(DATES == policy_dates{i}),':k',...
                  ['Testing phase: ' num2str(i)],...
                  'LabelVerticalAlignment','top',...
                  'LabelHorizontalAlignment','center',...
                  'LabelOrientation','aligned', ...
                  'interpreter','latex','HandleVisibility','off');
        end
    end


    %ylim([0 3e3])
    legend({'PHA 7-day' 'Model est.'},...
           'Location',[0.7 0.8 0.15 0.05],...
           'interpreter','latex')
    %title(['Symptomatic incidence'])
    hold off


    %title(['(R)ecovered' ' | ' region])
    ylabel('Incidence','interpreter','latex')
    ax = gca;
    ax.TickLabelInterpreter = 'latex';

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

    yyaxis left
    ylimits = ylim;
    ylim([0, min(max(datastore.(region).high95), 2*max(ydata))]);


    grid on
    box on

    yyaxis right
    tspan_tests=find(ismember(DATES,data_test.date));
    plot(t(tspan_tests),data_test.Nper100k(:,ri),'--','LineWidth',2,'HandleVisibility','off');
    %ylim([0 6])
    ylabel('Number of tests per 100 000','interpreter','latex')%'FontSize',8)
    hold off



    % save to file
    savepath = mfilename('fullpath');
    savepath1 = [savepath(1:end-15) 'Iinc_data_' region];
    if savetofile % true if we should compute the error table


        % finalize the plot
        h = gcf;
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        % polish target output size
        set(h,'PaperPositionMode','auto');
        set(h,'Position',[100 100 500 350]);
        print('-dpdf', savepath1)
        disp(['saved figure: ' savepath1]);
    else
        disp(['didn''t save figure: ' savepath1]);
    end

    %% Factor of underestimation (Model / Confirmed test)
    figure(2), clf

    % mu is the mean of the normal
    mu = exp(datastore.(region).muhat)';
    % sig is the std which can be found as the difference between:
    % sig = (mean) - (mean - 1 std)
    sig = exp(datastore.(region).muhat').*(1 - exp(-datastore.(region).sigmahat'));

    % Z = X / Y ~ Cauchy, but approx Normal with:
    mu_Z = mu ./ ydata(tspan_filter);
    delta_Y = ystd(tspan_filter) ./ ydata(tspan_filter); % CoV
    rho_Z = ystd(tspan_filter) ./ sig;
    sig_Z = sqrt(delta_Y.^2 .*( rho_Z.^(-2) + mu_Z));

    mu_Z = mu_Z';
    sig_Z = sig_Z';

    plot(t(tspan_filter),mu_Z,'Color', filter_c), hold on

    Ifac_68 = [mu_Z+sig_Z fliplr(mu_Z-sig_Z)];
    Ifac_95 = [mu_Z+2*sig_Z fliplr(mu_Z-2*sig_Z)];

    patch(pred_tt, Ifac_68,...
          [1 1 1],'FaceAlpha',0.15, ...
          'LineStyle','none','FaceColor',filter_c);
    patch(pred_tt, Ifac_95,...
          [1 1 1],'FaceAlpha',0.05, ...
          'LineStyle','none','FaceColor',filter_c);

    yline(1)

    % Testing policy change
    for i = 1:numel(policy_dates)
        if policy_dates{i} < DATES(end)
            xline(find(DATES == policy_dates{i}),':k',...
                  ['Testing phase: ' num2str(i)],...
                  'LabelVerticalAlignment','top',...
                  'LabelHorizontalAlignment','center',...
                  'LabelOrientation','aligned', ...
                  'interpreter','latex','HandleVisibility','off');
        end
    end



    %title(['(R)ecovered' ' | ' region])
    xtk = fliplr(tspan_filter(end:-28:1));
    %ylabel('Recovered [\%]','interpreter','latex')
    ax = gca;
    ax.TickLabelInterpreter = 'latex';

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

    ylabel('Est. symptomatic / Confirmed cases','interpreter','latex')

    % 2021 delimiter
    xline(find(DATES(tspan_filter) == 210101)+tspan_filter(1)-1,'--k','2021', ...
          'LabelVerticalAlignment','top','LabelOrientation','horizontal', ...
          'interpreter','latex','HandleVisibility','off');

    ylim([0 max(Ifac_68)])

    grid on
    box on
    hold off





    % save to file


    savepath = mfilename('fullpath');
    savepath2 = [savepath(1:end-15) 'Iinc_fac_' region];
    if savetofile % true if we should compute the error table


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


end



function ci = l_getci(lo95,lo68,mid,hi68,hi95,type)
    %L_GETCI formats the credible interval into the string pattern that we
    %  employ in the weekly report.
    %
    % lo95 - lower 95 limit
    % lo68 - lower 68 limit
    % mid  - mean
    % hi68 - upper 68 limit
    % hi95 - upper 95 limit
    % type - {'68'} or '95'

    ci = cell(size(lo68));
    switch type
      case '68'
        for i = 1:numel(lo68)
            ci{i} = ['[' num2str(lo68(i)) ','  num2str(mid(i)), ',' num2str(hi68(i)) ']'];
        end
      case '95'
        ci = ['[' num2str(lo95) ','  num2str(mid), ',' num2str(hi95) ']'];
      otherwise
        error('doesnt recognize the interval type')
    end
end

function tex = l_texify(TAB)
    %L_TEXIFY adds necessary additions into table TAB making it a bare minimum
    %   table in TeX.

    tex = '';
    tex = [tex '\begin{table}[htp]' newline ...
           '\centering' newline ...
           '\begin{tabular}{l ' repmat('r ',1,size(TAB,2)-1) '}' newline ...
           '\hline'];
    % per row
    for row = 1:size(TAB,1)
        ROW = '';
        % per column, append a '&' sign
        for col = 1:size(TAB,2)
            if col == size(TAB,2) % final column, add '\\' instead of '&'
                ROW = [ROW TAB{row,col} ' \\'];
            else
                ROW = [ROW TAB{row,col} '& '];
            end

        end % end cold
        tex = [tex newline ROW];
        if row == 1 % after first row, add additional line
            tex = [tex '\hline'];
        end
    end % end row

    tex = [tex newline ...
           '\hline' newline ...
           '\end{tabular}' newline ...
           '\caption{' ...
           'Uppskattad andel individer (\%) per region som genomgått ' ...
           'infektion (68\% KI samt medelvärde).}' ...
           '\textit{Notera:} för dessa estimat anses andelen individer ' ...
           'med genomgången infektion vara noll identiskt för startdatum ' ...
           '200401.' newline ...
           '\label{tab:recovered}' newline ...
           '\end{table}'];
end

function date = l_week2date(d_week,d_year)
    %L_WEEK2DATE converts a week number into a date, the Sunday of that week.

    t = datetime(d_year, 1, 10, 'Format', 'dd-MMMM-yyyy');  % First day of the year
    tspan = t - days(weekday(t)-1) + days((d_week-1).*7);

    yy = year(tspan);
    mm = month(tspan);
    dd = day(tspan);
    date = (yy-2000)*10000+mm*100+dd;

end
