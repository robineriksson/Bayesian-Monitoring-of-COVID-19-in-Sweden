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

% construct figures?
illustrate = true;

% which regions to include (see regionList below)
reg = [1];%1:21;

% dates of interest | for table
doi = [200801 201001 201201 210201 210401];
% ****************************************************


regionList = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
    'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
    'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
    'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
    'Jämtland' 'Västerbotten' 'Norrbotten'};

% matlab structs, does not allow åäö and ' ' in naming of fields.
regionList_nonnordic = strrep(regionList,'ö','o');
regionList_nonnordic = strrep(regionList_nonnordic,'Ö','O');
regionList_nonnordic = strrep(regionList_nonnordic,'ä','a');
regionList_nonnordic = strrep(regionList_nonnordic,'å','a');
regionList_nonnordic = strrep(regionList_nonnordic,' ','_');


datastore = struct();
datastore.table = [];
datastore.table_ci = [];
datastore.tablerows= {};

for i = reg
    region = regionList{i}; % to character
    region_ = regionList_nonnordic{i};


    % construct posterior file name
    % folder for posteriors
    abspath = mfilename('fullpath');
    prefix = [abspath(1:end-35) 'inference/results/SLAM/'];
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
    try
        load([abspath(1:end-35) 'weekly/save/runs/' posterior]);
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
    ixR = 7; % see weekly_prediction.m
    ZR = Z(ixR,:,:);
    stdZR = covZ.stdZ(ixR,:,:);



    %% Compute (R)ecovered.

    % Population data
    load Ncounties
    Npop = sum(N,1);
    regionList = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
        'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
        'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
        'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
        'Jämtland' 'Västerbotten' 'Norrbotten'};
    regi = find(ismember(regionList,region));

    YR = sum(ZR,1);
    YR = max(YR,1e-2);

    if strcmp(region,'Sweden')
        YR = YR / sum(Npop);
        stdYR = sum(stdZR,1) / sum(Npop);
    else
        YR = YR / Npop(regi);
        stdYR = stdZR / Npop(regi);
    end

    % Uncertainty
    sigmahat_R = sqrt(mean(log1p((stdYR./YR).^2),3));
    muhat_R = mean(log(YR),3);

    sigmahat_tot = sqrt(sigmahat_R.^2+var(log(YR),0,3));


    datastore.(region_).YR = YR;
    datastore.(region_).muhat = muhat_R;
    datastore.(region_).sigmahat = sigmahat_tot;
    datastore.(region_).DATES = DATES;
    datastore.(region_).tspan_filter = tspan_filter;

    %% process for table
    low95 = exp(datastore.(region_).muhat - 2*datastore.(region_).sigmahat);
    low68 = exp(datastore.(region_).muhat - 1*datastore.(region_).sigmahat);
    mid = exp(datastore.(region_).muhat);
    high68 = exp(datastore.(region_).muhat + 1*datastore.(region_).sigmahat);
    high95 = exp(datastore.(region_).muhat + 2*datastore.(region_).sigmahat);



    dates = DATES(tspan_data);
    t = 1:numel(DATES);

    % for abbreviations, see https://www.bydewey.com/monthdayabb.html
    if isempty(doi)
        strmonths = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'June'...
            'July' 'Aug' 'Sept' 'Oct' 'Nov' 'Dec'};
        % 1st each month:
        ixdates = find(mod(dates,100) == 1);
        % and final date
        ixdates = [ixdates; numel(dates)]; %+ final date
        xxdates = t(ixdates);
        % what month:
        months = mod(floor(dates/1e2),1e2);
        coltitle = strcat(strmonths(months(xxdates(1:end-1))), ' 1');

        coltitle = [coltitle, strcat(strmonths(months(xxdates(end))), ...
            [' ' num2str(mod(dates(end),100))])];
    else
        ixdates = find(ismember(dates,doi));
        % and final date
        ixdates = [ixdates; numel(dates)]; %+ final date
        xxdates = t(ixdates);
        coltitle = cellstr(num2str(dates(ixdates)))';
    end

    low95  = low95(ixdates);
    low68  = low68(ixdates);
    mid    = mid(ixdates);
    high68 = high68(ixdates);
    high95 = high95(ixdates);



    % store for the final table.
    datastore.table = cat(1, datastore.table, mid);
    % store the CI as well, for future reference.
    table_ci = sprintf('\\CI{%4.2f}{%4.2f}{%4.2f}{%4.2f}{%4.2f}\n',...
        [low95;low68;mid;high68;high95]*1e2);

    table_ci(end) = '';
    table_ci = split(table_ci,newline);
    datastore.table_ci = [datastore.table_ci, table_ci];


    %l_getci(low95,low68,mid,high68,high95,'68'));
    datastore.tablecols = coltitle;
    datastore.tablerows = cat(2,datastore.tablerows,region);
end
%%
% construct table (in latex format)
table_R = arr2latex(datastore.table,{'.2f'}, ...
    'collabel',{' ' datastore.tablecols{:}}, ...
    'rowlabel',datastore.tablerows,...
    'colspec',['l ' repmat('r ',1,numel(datastore.tablecols))],...
    'hline','off','centering','on', ...
    'caption','mean proportion of Recovered individuals', ...
    'label','tab:recovered');
%%
% add information to the _ci table
TAB = datastore.table_ci';
TAB = [datastore.tablecols; TAB];
TAB = [[{' '} regionList(reg)]', TAB];


TAB_tex = l_texify(TAB);

%%
% save .tex tables
tabname = [abspath(1:end-21) 'tab/tableRec.tex'];
if savetofile
    fileID = fopen(tabname,'w');
    fprintf(fileID,'%s\n',TAB_tex);
    fclose(fileID);
    disp(['saved table: ' tabname]);
else
    disp(['didn''t save table: ' tabname]);
end


%%

if illustrate

    %% load data from paper 1
    % Seropositivity in blood donors and pregnant women during the first year
    % of SARS-CoV-2 transmission in Stockholm, Sweden
    % https://onlinelibrary.wiley.com/doi/10.1111/joim.13304
    paper=csvread('RecPaperEstimate.csv',1);
    paper_dates=l_week2date(paper(:,1),paper(:,2));
    paper_pred = [paper_dates, paper(:,3:end)];
    xx_paper = find(ismember(dates,paper_dates));

        % at what dates do they mark measurements?
    measurements = [ (17:25),(30:34),(45:50),(4:8)
                    repmat(2020,1,20),repmat(2021,1,5)];
    paper_dates_m = l_week2date(measurements(1,:)',measurements(2,:)');
    xx_paper_m = find(ismember(dates,paper_dates_m));
    xx_paper_m_points = find(ismember(paper_dates,paper_dates_m));

    %% load data from study 1
    % Pavisning av antikroppar mot SARS-CoV-2 i blodprov från oppenvarden
    % https://www.folkhalsomyndigheten.se/contentassets/9c5893f84bd049e691562b9eeb0ca280/pavisning-antikroppar-mot-sars-cov-2-blodprov-oppenvarden.pdf
    %%
    study1=csvread('fhmAnti.csv',1);
    study1_dates=l_week2date(study1(:,2),study1(:,1));
    study1_pred= [study1_dates, study1(:,3:end)/100];
    % exlude final one, includes many vaccinated ind.
    study1_pred = study1_pred(1:end-1,:);
    xx_study1 = find(ismember(dates,study1_dates));

    %% load data from study 2
    % Påvisning av antikroppar mot SARS-CoV 2 hos blodgivare
    % https://www.folkhalsomyndigheten.se/contentassets/376f9021a4c84da08de18ac597284f0c/pavisning-antikroppar-mot-sars-cov-2-blodgivare.pdf
    study2=csvread('fhmAntiGivare.csv',1);
    study2_dates=l_week2date(study2(:,2),study2(:,1));
    study2_pred= [study2_dates, study2(:,3:end)/100];
    % exlude final one, includes many vaccinated ind.
    study2_pred = study2_pred(1:end-1,:);
    xx_study2 = find(ismember(dates,study2_dates));

    %% including vaccination numbers
    data_vac=loadData('FHM_vac');
    xx_vac = find(ismember(dates,data_vac.date));
    xx_vac_rev = find(ismember(data_vac.date,dates));
    %% Plot
    for i = reg
        region = regionList_nonnordic{i}
        figure(1), clf;
        % Sample trajectories
        Nsamples = 0;

          yy =squeeze(datastore.(region).YR(1,:,1:Nsamples));
          tspan_filter_short = datastore.(region).tspan_filter(1:end);

          % linear fit fit get a good start on 5 April
          lm = polyfit(1:10,paper_pred(1:10,4),1);
          start = polyval(lm,-0.7); % found this to work good.

          %
          yy = yy + start;
        if Nsamples > 0
          plot(tspan_filter_short,yy,...
            'Color',[0 0 1 0.1], 'HandleVisibility','off')
        end
        hold on


        % Mean
        mu = datastore.(region).muhat(1:end);
        sigma = datastore.(region).sigmahat(1:end);
        plot(tspan_filter_short,exp(mu)+start,'b')



        tt = [tspan_filter_short,fliplr(tspan_filter_short)];

        Y_R70 = [exp(mu-1.0364*sigma),fliplr(exp(mu+1.0364*sigma))] + start;
        Y_R95 = [exp(mu-1.9600*sigma),fliplr(exp(mu+1.9600*sigma))] + start;


        patch(tt, Y_R70,...
            [0.9 0.9 0.9],'FaceAlpha',0.15, ...
            'LineStyle','none','FaceColor','b');

        patch(tt,Y_R95, ...
            [0.9 0.9 0.9],'FaceAlpha',0.05, ...
            'LineStyle','none','FaceColor','b');

        % plot presented paper prediction
        plot(datastore.(region).tspan_filter(xx_paper),paper_pred(:,4),'r')
        plot(datastore.(region).tspan_filter(xx_paper_m),paper_pred(xx_paper_m_points,4),'or')

        pred_R70 = [paper_pred(:,3)', fliplr(paper_pred(:,5)')];
        pred_R95 = [paper_pred(:,2)', fliplr(paper_pred(:,6)')];
        pred_tt =  [datastore.(region).tspan_filter(xx_paper) ...
            fliplr(datastore.(region).tspan_filter(xx_paper))];
        patch(pred_tt, pred_R70,...
            [0.9 0.9 0.9],'FaceAlpha',0.15, ...
            'LineStyle','none','FaceColor','r');

        patch(pred_tt,pred_R95, ...
            [0.9 0.9 0.9],'FaceAlpha',0.05, ...
            'LineStyle','none','FaceColor','r');



        fhmcolor1 = [256, 0, 256];%[0 0 0];%
        fhmcolor2 = [100, 0, 256];%[0 0 0];%
        fhmalpha = 0.5;
        % plot presented FHM serology study 1
        plot(datastore.(region).tspan_filter(xx_study1),study1_pred(:,2),...
          'o','Color',[fhmcolor1 1]/256)
        for k = 1:size(xx_study1)
            x = datastore.(region).tspan_filter(xx_study1(k));
            linestyle='-';
            line([x x], [study1_pred(k,3), study1_pred(k,4)],...
              'LineStyle',linestyle,'Color',[fhmcolor1 fhmalpha*256]/256,'lineWidth',1.5)
        end

        % plot presented FHM serology study 2
        plot(datastore.(region).tspan_filter(xx_study2)+2,study2_pred(:,2),...
          'o','Color',[fhmcolor2 1]/256)
        for k = 1:size(xx_study2)
            x = datastore.(region).tspan_filter(xx_study2(k));
            linestyle='--';
            line([x+2 x+2], [study2_pred(k,3), study2_pred(k,4)],...
              'LineStyle',linestyle,'Color',[fhmcolor2 fhmalpha*256]/256,'lineWidth',1.5)
        end

%         study_R95 = [study_pred(:,3)', fliplr(study_pred(:,4)')];
%         pred_tt =  [datastore.(region).tspan_filter(xx_study) ...
%             fliplr(datastore.(region).tspan_filter(xx_study))];
%
%         patch(pred_tt,study_R95, ...
%             [0.9 0.9 0.9],'FaceAlpha',0.05, ...
%             'LineStyle','none','FaceColor','k');


        plot(datastore.(region).tspan_filter(xx_vac),data_vac.V1(xx_vac_rev,i)/Npop(i),':k')
        plot(datastore.(region).tspan_filter(xx_vac),data_vac.V2(xx_vac_rev,i)/Npop(i),'-.k')


        pos = 5;
        ylims = [0 0.3];


        a1 = annotation('textarrow',...
          [datastore.(region).tspan_filter(xx_vac(pos))*0.82 datastore.(region).tspan_filter(xx_vac(pos))*0.90]/tspan_data(end),...
          [data_vac.V1(xx_vac_rev(pos),i)/Npop(i)*3.8, data_vac.V1(xx_vac_rev(pos),i)/Npop(i)*3]/ylims(end),...
          'String','\%1st dose', 'interpreter','latex');
        a1.LineStyle=':';
        a2 = annotation('textarrow',...
          1.1*[datastore.(region).tspan_filter(xx_vac(pos))*0.96 datastore.(region).tspan_filter(xx_vac(pos))*0.88]/tspan_data(end),...
          15*[data_vac.V2(xx_vac_rev(pos),i)/Npop(i)*3.3, data_vac.V2(xx_vac_rev(pos),i)/Npop(i)*3]/ylims(end),...
          'String','\%2nd dose', 'interpreter','latex');
        a2.LineStyle='-.';

        hold off

        axis([tspan_data([1 end]) ylims]);

        %title(['(R)ecovered' ' | ' region])
        xtk = fliplr(tspan_data(end:-28:1));
        ylabel('Recovered [\%]','interpreter','latex')
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
%         set(gca,'xtick',xtk);
%
%         xlabs = string(DATES(xtk));
%         % thin the xticks
%         xlabs(2:2:end) = "";
%
%
%
%         set(gca,'xticklabel',xlabs);
%         xtickangle(45);
%

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



        grid on
        box on
        leg = {'Median','70\%','95\%'};%,'[JIM 2021] Median','70\%','95\%'};
        legend(leg, ...
            'Location','northwest','interpreter','latex','FontSize',12,...
            'Orientation','horizontal',...
            'NumColumns', 1);

        % change y to \%
        yticklabels(num2cell(100*yticks))
        %%
        % save
        figname = [abspath(1:end-21) 'fig/Recov_' region];

        if savetofile
            h = gcf;

            set(h,'PaperPositionMode','auto');
            set(h,'Position',[100 100 500 350]);

            print('-dpdf', figname)

            disp(['saved figure: ' figname])
        else
          disp(['didn''t save figure: ' figname])
        end

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
