% ARX_FIT constructs an AR prediction model using an expanding window
% formulation. The idea of the script is to generate a model for
% comparison with the posterior Kalman filter. How does the prediction
% differ and so on. The comparison between the models is done by
% loading the results from previous weekly reportings with
% READWEEKLYTABLE. We can therefore evaluate the prediciton inside 66/95% CrI
% and NRMSE) of the AR and the Posterior Kalman filter on the same
% dates, thus making sure not to favorize the Posterior Kalman filter
% by "evaluating on already seen data".
%
% The AR model uses the ARX function supplied by MATLAB'S libraries.
%
% The AR model is basic by design. Further effort could have been
% included, by for example extending the model with the the ARMA or
% ARIMA extensions.

% R. Marin 2022-04-07
if ~exist('savetofile','var')
    savetofile=false;
end
if ~exist('reg','var')
    reg=2;
end
if ~exist('verb','var')
    verb=false;
end

if ~exist('alldata','var')
    alldate=true;
end

if ~exist('evalconst','var')
    evalconst=false;
end
inferDincPlot = false; % true if generate Dinc as diff(D)

register = 'C19'; % data source

startdate = 200401;
enddate = 210531;
kalmancompare = 0;
% load and clean data
data_all = loadData(register);
indata = find(data_all.date >= startdate & data_all.date <= enddate);

start=1;
data_all.D = data_all.D(indata,:);
data_all.Dinc = data_all.Dinc(indata,:);
data_all.H = data_all.H(indata,:);
data_all.Icum = data_all.Icum(indata,:);
data_all.W = data_all.W(indata,:);
data_all.date = data_all.date(indata);

data_all = polishData(data_all,'D','Dinc',1);
data_all = smoothData(data_all,{'D' 'H' 'W'},{'Dinc' [] []});

if ~exist('ypred','var')


    len = size(data_all.date,1);


    stopvec = 10:len;
    ypred= nan(len,3);
    ysdpred = nan(len,3);

    if ~alldata
        stopvec = [251 258 265 272 279 286 293 307 314 321 ...
                   328 335 342 349 356 363 370 377 384 391 ...
                   398 405 412 419 426];
    end


    % expanding window formulation
    if verb
        disp('fitting AR ...');
    end
    for stop = stopvec
        if mod(stop,20)==0
            if verb
                disp([num2str(round(100*stop/len)) '%'])
            end
        end
        z = iddata([data_all.H(1:stop,reg), data_all.W(1:stop,reg),...
                    data_all.D(1:stop,reg)],[],1,'TimeUnit','days','Tstart',1);


        % delinearize for better fit
        [zd, Tz] = detrend(z, 0);

        % constans choosen for ~good predictions the first 50 days
        na = [1 1 1]';
        nb = [0 1 0]';
        nc = [1 1 1]';

        % multi variance ARMA model
        sysARMA = arx(zd,[na nb nc]);

        % do a seven day ahead prediction
        forecastOpt = forecastOptions('OutputOffset',Tz.OutputOffset');
        [yf1,~,~,ysd1] = forecast(sysARMA, z, 7, forecastOpt);

        % and save the 7th day
        ypred(stop,:) = yf1.y(end,:);
        ysdpred(stop,:) = ysd1(end,:);
    end
    ydata = z.y;
end



%% evaluate
TAB68 = zeros(2,3);
TAB95 = zeros(2,3);

NRMSE = @(x,y) sqrt(nanmean( (x-y).^2,2))./nanmean(x,2);

% Selected Mondays predictions
% -load Kalman predictions
[kalman_reported, date_comp, TABkalman] = l_getKalmanScore();
% kalman_reported = {{data} {prediction}}

TAB68(1,:) = 1-TABkalman(1,:,1);
TAB95(1,:) = 1-TABkalman(1,:,2);

% nrmse for Kalman
kalman_nrmse = NRMSE(kalman_reported{1}.reported,kalman_reported{2}.reported)';

% unpack for energyscore
kdata = [kalman_reported{1}.H',kalman_reported{1}.W',kalman_reported{1}.D'];
kpred = [kalman_reported{2}.H',kalman_reported{2}.W',kalman_reported{2}.D'];
ksdpred = [kalman_reported{2}.Hsd',kalman_reported{2}.Wsd',kalman_reported{2}.Dsd'];
kalman_es = energyscore(kdata,kpred,ksdpred,1e4);




t_reported = find(ismember(data_all.date(:,end),date_comp));

% inside | outside?
Yhi95 = ypred + 1.96*ysdpred;
Ylo95 = ypred - 1.96*ysdpred;
Yhi68 = ypred + 1.00*ysdpred;
Ylo68 = ypred - 1.00*ysdpred;

% find outside
elem_reported95 = [Yhi95(t_reported,1)' < ydata(t_reported,1)' | Ylo95(t_reported,1)' > ydata(t_reported,1)',
                   Yhi95(t_reported,2)' < ydata(t_reported,2)' | Ylo95(t_reported,2)' > ydata(t_reported,2)',
                   Yhi95(t_reported,3)' < ydata(t_reported,3)' | Ylo95(t_reported,3)' > ydata(t_reported,3)'];
elem_reported68 = [Yhi68(t_reported,1)' < ydata(t_reported,1)' | Ylo68(t_reported,1)' > ydata(t_reported,1)',
                   Yhi68(t_reported,2)' < ydata(t_reported,2)' | Ylo68(t_reported,2)' > ydata(t_reported,2)',
                   Yhi68(t_reported,3)' < ydata(t_reported,3)' | Ylo68(t_reported,3)' > ydata(t_reported,3)'];

% report inside  = 1-outside
TAB68(2,:) = 1-mean(elem_reported68,2);
TAB95(2,:) = 1-mean(elem_reported95,2);


% compute the (mean normalized) NRMSE
arx_nrmse = NRMSE(ypred(t_reported,:)',ydata(t_reported,:)')';

% and the Energy score
arx_es = energyscore(ydata(t_reported,:),ypred(t_reported,:),ysdpred(t_reported,:),1e4);

TAB68 = cat(1,[{'Hospital (H)'} {'Intensive (W)'} {'Death (D)'}],num2cell(TAB68));
TAB68 = cat(2,[{' '} {'Posterior Kalman (68\% CrI)'} {'AR'}]',TAB68);
TAB95 = cat(2,[{'Posterior Kalman (95\% CrI)'} {'AR'}]',num2cell(TAB95));
TAB = cat(1,TAB68,TAB95);

TAB_nrmse = cat(2,[{'Posterior Kalman (NRMSE)'} {'AR'}]',num2cell([kalman_nrmse;arx_nrmse]));
TAB = cat(1,TAB,TAB_nrmse);

TAB_es = cat(2,[{'Posterior Kalman (Energy score)'} {'AR'}]',num2cell([kalman_es;arx_es]./100));
if verb
    TAB = cat(1,TAB,TAB_es)
end


tablePoly = l_latexify(TAB,numel(t_reported),[3 5 7]);
if verb
    tablePoly
end

%% plotting
ypred_plot = ypred;
ysdpred_plot = ysdpred;
ydata_plot = ydata;


smoothing = @(x) sgolayfilt(x,0,7);

if inferDincPlot
    ydata_plot(:,3) = smoothing([0;diff(ydata_plot(:,3),1)]);
else
    ydata_plot(:,3) = smoothing(data_all.Dinc(:,2));
end

ypred_plot(:,3) = smoothing([0;diff(ypred_plot(:,3),1)]);
ysdpred_plot(:,3) = [0; sqrt(ysdpred(2:end,3).^2+ysdpred(1:end-1,3).^2)];

color = [[0 0 1]; % H
         [1 0 0]; % W
         [160 5 30]/256]; % D

figure(1), clf
yyaxis left, hold on
plot(ydata_plot(:,1),'.','Color',color(1,:),'HandleVisibility','off')
yyaxis right, hold on
plot(ydata_plot(:,2),'.','Color',color(2,:),'HandleVisibility','off')
plot(ydata_plot(:,3),'.','Color',color(3,:),'HandleVisibility','off')

tt = [[1:len] fliplr(1:len)];
yy68 = [(ypred_plot + ysdpred_plot)' fliplr((ypred_plot - ysdpred_plot)')];



for i = 1:3
    if ismember(i,[2 3])
        yyaxis right
    else
        yyaxis left
    end
    notnan = ~isnan(yy68(i,:));
    patch(tt(notnan),yy68(i,notnan), ...
          [0.9 0.9 0.9],'FaceAlpha',0.15, ...
          'LineStyle','None', 'FaceColor',color(i,:),...
          'HandleVisibility','off');
    plot(ypred_plot(:,i),'-','Color',color(i,:))
end

xline(t_reported(1),'HandleVisibility','off','Label','Reporting',...
     'interpreter','latex');
patch([t_reported(1) t_reported(end) t_reported(end) t_reported(1)],...
      [0 0 32 32],[0.9 0.9 0.9],...
      'FaceAlpha',0.075, ...
      'LineStyle','None', 'FaceColor',[0 0 0],...
      'HandleVisibility','off');



axis tight
grid on
yyaxis left
ylim([0 160])
yyaxis right
ylim([0 32])
yticks(0:4:32)

% update ticks
TSPAN = z.SamplingInstants;
DATES = datenum(2e3+floor(data_all.date(1)/1e4), ...
                mod(floor(data_all.date(1)/1e2),1e2) ,...
                mod(data_all.date(1),1e2));
DATES = DATES:DATES+numel(TSPAN)-1;
DATES = datevec(DATES);
DATES = DATES(:,1:3)*[1e4 1e2 1]'-2e7; % finally!
                                       %tspan_data = TSPAN(x_pred);  % data used
tspan_data = TSPAN;  % data used



xtk = TSPAN(start+5:28:end);
grid on;
box on;

% wanted xticks
% for abbreviations, see https://www.bydewey.com/monthdayabb.html
strmonths = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'June'...
             'July' 'Aug' 'Sept' 'Oct' 'Nov' 'Dec'};

% 15th each month:
ixdates = find(mod(DATES(tspan_data),100) == 1);
xxdates = tspan_data(ixdates);
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
xline(find(DATES(tspan_data) == 210101)+tspan_data(1)-1,'--k','2021', ...
      'LabelVerticalAlignment','top','LabelOrientation','horizontal', ...
      'interpreter','latex','HandleVisibility','off');

% legends
leg = {'Hospital [H]','Intensive [W] $\,$(right)',...
       'Dead/day [D] (right)'};

legend(leg, ...
       'Location','northwest','interpreter','latex','FontSize',12,...
       'Orientation','horizontal',...
       'NumColumns', 1);



h = gcf;
ax = gca;
ax.TickLabelInterpreter = 'latex';


hold off


savepath = mfilename('fullpath');
savepath = savepath(1:end-15);

tabname = [savepath '../tab/' 'tablePolyfit.tex'];
figname1 = [savepath 'arx.pdf'];

if savetofile
  fileID = fopen(tabname,'w');
  fprintf(fileID,'%s\n',tablePoly);
  fclose(fileID);
  if verb
      disp(['saved table: ' tabname]);
  end
  figure(1)
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 500 350]);
  print(gcf, figname1, '-dpdf')
  if verb
      disp(['saved figure: ' figname1]);
  end
else
    if verb
        disp(['didn''t save table: ' tabname]);
        disp(['didn''t save figure: ' figname1]);
    end
end


if evalconst
    figure(2)
    evalcorr(zd.y,50);
end

function [kalman_reported, date_comp, TABkalman] = l_getKalmanScore()
%L_GETKALMANSCORE unpacks the READWEEKLYTABLE output into the format we want
    [Uppsala,~,~,week,date_comp,TABkalman] = readWeeklyTable();

    % Compared to NRMSE on Kalman on the selected predictions
    % - but first we need to extract the predictions and data.
    kalman_pred = struct('H',zeros(1,numel(Uppsala)),...
                         'W',zeros(1,numel(Uppsala)),...
                         'D',zeros(1,numel(Uppsala)),...
                         'Hsd',zeros(1,numel(Uppsala)),...
                         'Wsd',zeros(1,numel(Uppsala)),...
                         'Dsd',zeros(1,numel(Uppsala)));

    kalman_data = struct('H',zeros(1,numel(Uppsala)),...
                         'W',zeros(1,numel(Uppsala)),...
                         'D',zeros(1,numel(Uppsala)));

    for k = 1:numel(Uppsala)
        % mean predictions
        kalman_pred.H(k) = Uppsala{k}(1,3);
        kalman_pred.W(k) = Uppsala{k}(2,3);
        kalman_pred.D(k) = Uppsala{k}(3,3);

        kalman_pred.Hsd(k) = mean([Uppsala{k}(1,3)-Uppsala{k}(1,2), Uppsala{k}(1,4)-Uppsala{k}(1,3)]);
        kalman_pred.Wsd(k) = mean([Uppsala{k}(2,3)-Uppsala{k}(2,2), Uppsala{k}(2,4)-Uppsala{k}(2,3)]);
        kalman_pred.Dsd(k) = mean([Uppsala{k}(3,3)-Uppsala{k}(3,2), Uppsala{k}(3,4)-Uppsala{k}(3,3)]);

        % the data at that time | sadly data has been updated cont.
        kalman_data.H(k) = Uppsala{k}(1,end);
        kalman_data.W(k) = Uppsala{k}(2,end);
        kalman_data.D(k) = Uppsala{k}(3,end);
    end

    kalman_data.reported = [kalman_data.H;kalman_data.W;kalman_data.D];
    kalman_pred.reported = [kalman_pred.H;kalman_pred.W;kalman_pred.D];

    kalman_reported = {kalman_data, kalman_pred};
end

function tex = l_latexify(TAB,N,hmid)
%L_LATEXIFY generates a table environment with a tabular of the
% correct format for the paper result.


    X = cell2mat(TAB(2:end,2:end));

    %caption = ['RMSE of weekly predictions for proposed model and baseline models'];
    caption = ['Frequency of all weekly reported predictions (N~=~' num2str(N) ','...
               ' Dec 2020--May 2021) for Uppsala that fell inside of the reported CrIs' ...
               ' (68/95\%, 7 days ahead), the NRMSE, and the multivariate Energy Score' ...
               ' evaluated on the following week,' ...
               ' thus comparing the performance of the posterior Kalman filter and the AR model.'];


    label = ['tab:modelcomparison'];
    colspec = ['l' 'r' 'r' 'r'];

    tex ='';
    tex = [tex '\begin{table}[htp]' newline ...
           '\centering' newline ...
           '\caption{' caption '}' newline ...
           '\label{' label '}' newline ...
           '\begin{tabular}{' colspec '}'];
    tex = [tex '\hline'];
    % per row
    for row = 1:size(TAB,1)
        ROW = '';
        % per column, append a '&' sign
        for col = 1:size(TAB,2)
            val = TAB{row,col};
            if isnumeric(val)
                val = num2str(100*round(val,2,'significant'));
            end



            if col == size(TAB,2) % final column, add '\\' instead of '&'
                ROW = [ROW val '\\'];
                if ismember(row,hmid)
                    ROW = [ROW '\hline'];
                end

            else
                ROW = [ROW val '& '];
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
           '\end{table}'];
end

function [] = evalcorr(y,stop)
% EVALCORR generates the acf and pacf plots plots per data dimensions
% to help evaluate the order of AR constants. y is the input data, and
% needs to be of size [N x 3]. STOP indicates the datastop for
% training.
%
% A systematic approach is to use the following programatic solutions
% given by Matlab:
% V1 = arxstruc(y(1:stop-10,1),zd(stop+9:stop,1),(1:10)');
% V2 = arxstruc(y(1:stop-10,2),zd(stop+9:stop,2),(1:10)');
% V3 = arxstruc(y(1:stop-10,3),zd(stop+9:stop,3),(1:10)');
% na1 = selstruc(V1,0)
% na2 = selstruc(V2,0)
% na3 = selstruc(V3,0)

    [n,m] = size(y);

    stop = min(n,stop);

    for i = 1:3
        subplot(2,3,i)
        autocorr(y(1:stop,i),stop-1)
    end

    for i = 4:6
        subplot(2,3,i)
        parcorr(y(1:stop,i-3),stop-1)
    end
end
