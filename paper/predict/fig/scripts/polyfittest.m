%POLYFITTEST runs a weighted linear model fit on the data (rolling 30 day
% window) and runs weekly predictions.
%

% R. Eriksson 12-29-2021

%% Prepare data
% control variables
lookback_window = 30; % N-days rolling window
pred_window = 7; % predict 1-M days ahead

% Note about kappa, this could (should) be choosen by an optimization
% schema that find the optimum between bias and variance. An overfit value
% could be found men setting kappa high ~ 0.5, it's jumpy and a less greedy
% approach is a small kappa ~ 0.01, which is similar to the non weighted
% approach.
kappa = 0.15; % weight scaling [w = exp(kappa*x)]

fitalpha95 = 0.05;% 95%; % linear model CI level [100(1-alpha)]
fitalpha68 = 0.32; % 68%
register = 'C19'; % data source
reg = 'Uppsala'; % data region
startdate = 200401;
enddate = 210531;
kalmancompare = 0;%load in full kalman data or not. (if RMSE on all dates is wanted)
% if exporting the table to latex and saving it together with the two
% figures.
save2file=false;



% load and clean data
data_all = loadData(register);
indata = find(data_all.date >= startdate & data_all.date <= enddate);
%indata = find(data_all.date <= enddate);
%start = find(data_all.date >= startdate, 1);
start=1;
data_all.D = data_all.D(indata,:);
data_all.Dinc = data_all.Dinc(indata,:);
data_all.H = data_all.H(indata,:);
data_all.Icum = data_all.Icum(indata,:);
data_all.W = data_all.W(indata,:);
data_all.date = data_all.date(indata);

data_all = polishData(data_all,'D','Dinc',1);
data_all = smoothData(data_all,{'D' 'H' 'W'},{'Dinc' [] []});

%% RMSE (and UC check) for the linear filter
nrmse_lm_all = zeros(3,0);
nrmse_lm_reported = zeros(3,0);
TAB95 = cell(3,3);
TAB68 = cell(3,3);

% error function
NRMSE = @(x,y) sqrt(mean( (x-y).^2,2))./mean(x,2);


% Selected Mondays predictions
% -load Kalman predictions
[Uppsala,~,~,week,date_comp,TABkalman] = readWeeklyTable();


% Compared to NRMSE on Kalman on the selected predictions
% - but first we need to extract the predictions and data.
kalman_pred = struct('H',zeros(1,numel(Uppsala)),...
  'W',zeros(1,numel(Uppsala)),...
  'D',zeros(1,numel(Uppsala)));
kalman_data = struct('H',zeros(1,numel(Uppsala)),...
  'W',zeros(1,numel(Uppsala)),...
  'D',zeros(1,numel(Uppsala)));

for k = 1:numel(Uppsala)
  % mean predictions
  kalman_pred.H(k) = Uppsala{k}(1,3);
  kalman_pred.W(k) = Uppsala{k}(2,3);
  kalman_pred.D(k) = Uppsala{k}(3,3);

  % the data at that time | sadly data has been updated cont.
  kalman_data.H(k) = Uppsala{k}(1,end);
  kalman_data.W(k) = Uppsala{k}(2,end);
  kalman_data.D(k) = Uppsala{k}(3,end);
end

kalman_data.reported = [kalman_data.H;kalman_data.W;kalman_data.D];
kalman_pred.reported = [kalman_pred.H;kalman_pred.W;kalman_pred.D];
% % flip the ordering for easier plotting
% fn = fieldnames(kalman_pred);
% for k=1:numel(fn)
%   kalman_pred.(fn{k}) = fliplr(kalman_pred.(fn{k}));
%   kalman_data.(fn{k}) = fliplr(kalman_data.(fn{k}));
% end


nrmse_kalman_reported = NRMSE(kalman_pred.reported,kalman_data.reported)


if kalmancompare
  lagplot % to get muhat and Ydata_swe
  %%
  nrmse_kalman_all = NRMSE(exp(muhat(:,1:end-7)),Ydata_swe);
else
  nrmse_kalman_all = NaN;
end

%TAB(1,:) = [{sum(rmse_kalman_reported)} {sum(rmse_kalman_all)} {NaN}];
TAB68(1,:) = [num2cell(1-TABkalman(1,:,1))];
TAB95(1,:) = [num2cell(1-TABkalman(1,:,2))];


%% Linear model (simple + exponential)
for useweights = [0 1] % 0: standard poly, 1: exp weighted poly
  % rolling lookback window
  regid = strcmp(reg,regions);
  data = struct('H', l_rolling(data_all.H(:,regid),lookback_window),...
    'W', l_rolling(data_all.W(:,regid),lookback_window),...
    'D', l_rolling(data_all.D(:,regid),lookback_window),...
    'date', l_rolling(data_all.date,lookback_window),...
    'x',l_rolling(1:numel(data_all.date),lookback_window));

  % from numeric to datetime style
  lag = 0; % 0 if use full prediction window in eval, else > 0
  %x_pred = data.x(:,end)+pred_window-lag;
  x_pred = max(data.x,[],2)-lag;
  x_data = 1:numel(data_all.date);

  TSPAN = x_data(1):x_pred(end);
  DATES = datenum(2e3+floor(data_all.date(1)/1e4), ...
    mod(floor(data_all.date(1)/1e2),1e2) ,...
    mod(data_all.date(1),1e2));
  DATES = DATES:DATES+numel(TSPAN)-1;
  DATES = datevec(DATES);
  DATES = DATES(:,1:3)*[1e4 1e2 1]'-2e7; % finally!
  %tspan_data = TSPAN(x_pred);  % data used
  tspan_data = TSPAN(x_pred);  % data used

  % prepare prediction data
  rows = size(data.H,1);
  val = struct('H', zeros(rows, size(data.H,2) + pred_window),...
    'Hci95', zeros(rows, size(data.H,2) + pred_window, 2),...
    'Hci68', zeros(rows, size(data.H,2) + pred_window, 2),...
    'W', zeros(rows, size(data.W,2) + pred_window),...
    'Wci95', zeros(rows, size(data.H,2) + pred_window,2),...
    'Wci68', zeros(rows, size(data.H,2) + pred_window,2),...
    'D', zeros(rows, size(data.D,2) + pred_window),...
    'Dci95', zeros(rows, size(data.H,2) + pred_window,2),...
    'Dci68', zeros(rows, size(data.H,2) + pred_window,2),...
    'x', zeros(rows, size(data.D,2) + pred_window));




  % exponentially decaying weights
  if useweights
    disp('using weights')
    weights = exp(kappa*(1:size(data.H,2))); % exponentially increasing weights
  else
    disp('not using weights')
    weights = ones(1,size(data.H,2));
  end
  weights = weights / sum(weights); % normalize weights

  %% fit and predict
  for k = 1:rows
      % fit
      lm_H=fitlm(tbl,'H~x+W+D')
      lm_W=fitlm(tbl,'W~x+H+D')
      lm_D=fitlm(tbl,'D~x+H+W')

      poly = struct('H', fitlm(data.x(k,:),data.H(k,:), 'Weights',weights),...
                    'W', fitlm(data.x(k,:),data.W(k,:), 'Weights',weights),...
                    'D', fitlm(data.x(k,:),data.D(k,:), 'Weights',weights));
      % predict
      val.x(k,:) = [data.x(k,:), (1:pred_window) + max(data.x(k,:))];

      [val.H(k,:), val.Hci95(k,:,:)] = predict(poly.H, val.x(k,:)', 'Alpha', fitalpha95,'Prediction','observation');
      [val.W(k,:), val.Wci95(k,:,:)] = predict(poly.W, val.x(k,:)', 'Alpha', fitalpha95,'Prediction','observation');
      [val.D(k,:), val.Dci95(k,:,:)] = predict(poly.D, val.x(k,:)', 'Alpha', fitalpha95,'Prediction','observation');

      [~, val.Hci68(k,:,:)] = predict(poly.H, val.x(k,:)', 'Alpha', fitalpha68,'Prediction','observation');
      [~, val.Wci68(k,:,:)] = predict(poly.W, val.x(k,:)', 'Alpha', fitalpha68,'Prediction','observation');
      [~, val.Dci68(k,:,:)] = predict(poly.D, val.x(k,:)', 'Alpha', fitalpha68,'Prediction','observation');

  end

  % Make sure predictions are > 0
  val.H(val.H < 0) = 0;
  val.Hci95(val.Hci95 < 0) = 0;
  val.Hci68(val.Hci68 < 0) = 0;
  val.W(val.W < 0) = 0;
  val.Wci95(val.Wci95 < 0) = 0;
  val.Wci68(val.Wci68 < 0) = 0;
  val.D(val.D < 0) = 0;
  val.Dci95(val.Dci95 < 0) = 0;
  val.Dci68(val.Dci68 < 0) = 0;


  %% visually evaluate the predictions
  Hcomp = 1; Wcomp = 2; Dcomp = 3; % bookkeeping
  color = [[0 0 1]; % H
    [1 0 0]; % W
    [160 5 30]/256; % D
    [0 1 0]]; % I


  %
  figure(useweights+1), clf;

  yyaxis left % Hospitalized
  plot(x_pred,val.H(:,end-lag),'-', 'Color', color(Hcomp,:)), hold on
  plot(x_data,data_all.H(:,regid),'.', 'Color', color(Hcomp,:), ...
    'HandleVisibility','off')

  if strcmp(reg,'Uppsala')
    ylim([0,160]) % match figure in paper
    yticks(0:20:160)
  end


  yyaxis right % Worse (ICU) & Dead (per day)
  plot(x_pred,val.W(:,end-lag),'-', 'Color', color(Wcomp,:))
  plot(x_data,data_all.W(:,regid), '.','Color', color(Wcomp,:), ...
    'HandleVisibility','off')


  smoothing = @(x) sgolayfilt(x,0,7);
  plot(x_pred,smoothing([0;diff(val.D(:,end-lag),1)]),'-', ...
    'Color', color(Dcomp,:))
  plot(x_data,smoothing([0; diff(data_all.D(:,regid))]), '.', ...
    'Color', color(Dcomp,:), 'HandleVisibility','off')

  if strcmp(reg,'Uppsala')
    ylim([0,32]) % same as above
    yticks(0:4:32)
  end

  % shaded region (confidence interval)
  tt = [x_pred' fliplr(x_pred')];
  yy95 = [squeeze(val.Hci95(:,end-lag,2))' fliplr(squeeze(val.Hci95(:,end-lag,1))')];
  yy95 = cat(1,yy95,...
    [squeeze(val.Wci95(:,end-lag,2))' fliplr(squeeze(val.Wci95(:,end-lag,1))')]);
  DincCI95 = smoothing([[0 0];diff(squeeze(val.Dci95(:,end-lag,:)),1)]);
  yy95 = cat(1,yy95,[DincCI95(:,1)' fliplr(DincCI95(:,2)')]);

  for i = 1:3
    % note: lognormal interpretation
    if any(i == [Hcomp])
      yyaxis left
    else
      yyaxis right
    end

    % when NaN values are included, it cannot fill, only supplies
    % lines. Repalce with zeros.
    [row, col] = find(isnan(yy95));
    yy95(row,col) = 0;
    patch(tt,yy95(i,:), ...
      [0.9 0.9 0.9],'FaceAlpha',0.15, ...
      'LineStyle','none', 'FaceColor',color(i,:),...
      'HandleVisibility','off');
  end

  %%

  % update ticks
  xlim([start,x_pred(end-7)])
  xtk = TSPAN(start+5:28:end);
  grid on;

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
  leg = {'Hospital [H]','Intensive [W]','Dead/day [D]'};

legend(leg, ...
  'Location','northwest','interpreter','latex','FontSize',12,...
  'Orientation','horizontal',...
  'NumColumns', 1);



  h = gcf;
  ax = gca;
  ax.TickLabelInterpreter = 'latex';


  hold off


  %% Eval the predictions

  % on all one week ahead predictions
  [row, col] = max(data.x(8:end,:),[],2);
  H=NaN(1,numel(row));
  W=NaN(1,numel(row));
  D=NaN(1,numel(row));
  %x=zeros(1,numel(row));

  for k = 1:numel(row)
      H(k) = data.H(k+7,col(k));
      W(k) = data.W(k+7,col(k));
      D(k) = data.D(k+7,col(k));
      %x(k) = data.x(k+7,col(k));
  end
  %sum(x' - val.x(1:end-7,end)) % checks the date match
  nrmse_lm_all = cat(2,nrmse_lm_all,...
                    [NRMSE(H, val.H(1:end-7,end)')
                     NRMSE(W, val.W(1:end-7,end)')
                     NRMSE(D, val.D(1:end-7,end)')]);


  % find the matching dates
  data_reported = find(ismember(data.date(:,end),date_comp));
  pred_reported = [data_reported(1)-7;  data_reported(1:end-1)];

  %sum(data.x(data_sel,end) - val.x(pred_sel,end)) % checks the date match
  nrmse_lm_reported = cat(2,nrmse_lm_reported,...
                    [NRMSE(data.H(data_reported,end)', val.H(pred_reported,end)')
                     NRMSE(data.W(data_reported,end)', val.W(pred_reported,end)')
                     NRMSE(data.D(data_reported,end)', val.D(pred_reported,end)')]);

  %TAB(useweights+2,:) = [{sum(rmse_lm_reported)} {sum(rmse_lm_all)}];

  % inside or outsize?
%   elem_all = [val.Hci95(1:end-7,end-lag,1)' > H | val.Hci95(1:end-7,end-lag,2)'<H
%               val.Wci95(1:end-7,end-lag,1)' > W | val.Wci95(1:end-7,end-lag,2)'<W
%               val.Dci95(1:end-7,end-lag,1)' > D | val.Dci95(1:end-7,end-lag,2)'<D];
%   inside_all = 1-mean(elem_all,2);

  elem_reported95 = [val.Hci95(pred_reported,end,1)' > data.H(data_reported,end)' | val.Hci95(pred_reported,end,2)' < data.H(data_reported,end)'
              val.Wci95(pred_reported,end,1)' > data.W(data_reported,end)' | val.Wci95(pred_reported,end,2)' < data.W(data_reported,end)'
              val.Dci95(pred_reported,end,1)' > data.D(data_reported,end)' | val.Dci95(pred_reported,end,2)' < data.D(data_reported,end)'];
  inside_reported95 = 1-mean(elem_reported95,2);

  elem_reported68 = [val.Hci68(pred_reported,end,1)' > data.H(data_reported,end)' | val.Hci68(pred_reported,end,2)' < data.H(data_reported,end)'
              val.Wci68(pred_reported,end,1)' > data.W(data_reported,end)' | val.Wci68(pred_reported,end,2)' < data.W(data_reported,end)'
              val.Dci68(pred_reported,end,1)' > data.D(data_reported,end)' | val.Dci68(pred_reported,end,2)' < data.D(data_reported,end)'];
  inside_reported68 = 1-mean(elem_reported68,2);

  TAB95(useweights+2,:) = num2cell(inside_reported95)';
  TAB68(useweights+2,:) = num2cell(inside_reported68)';

end
%%
% %%
% figure(2)
% subplot(1,3,1)
% plot(data_sel,[data.H(data_sel,end)'; kalman_data.H], '.'), hold on
% plot(pred_sel, val.H(pred_sel,end)','b')
% plot(data_sel,kalman_pred.H,'r'), hold off
% ylabel('H')
%
% set(gca,'xticklabel',string(fliplr(week)));
%
% subplot(1,3,2)
% plot(data_sel,[data.W(data_sel,end)'; kalman_data.W], '.'), hold on
% plot(pred_sel, val.W(pred_sel,end)','b')
% plot(data_sel,kalman_pred.W,'r'), hold off
% ylabel('W')
% xlabel('week')
% set(gca,'xticklabel',string(fliplr(week)));
%
% subplot(1,3,3)
% plot(data_sel,[data.D(data_sel,end)'; kalman_data.D], '.'), hold on
% plot(pred_sel, val.D(pred_sel,end)','b')
% plot(data_sel,kalman_pred.D,'r'), hold off
% ylabel('D')
% set(gca,'xticklabel',string(fliplr(week)));
% legend({'lm data', 'kalman data', 'lm pred', 'kalman pred'},...
% 'Location','southeast')

%%
nrmse_fit = [nrmse_kalman_reported, nrmse_lm_reported]
nrmse_fit_mean = sqrt(mean(nrmse_fit.^2,1))

%% Build the Table
%TAB = cat(1,[{'RMSE all'} {'RMSE reporting dates'}],TAB);
%TAB = cat(1,[{'Inside all'} {'Inside reporting dates'}],TAB);
TAB68 = cat(1,[{'Hospital (H)'} {'Intensive (W)'} {'Death (D)'}],TAB68);
TAB68 = cat(2,[{' '} {'Posterior Kalman (68\%)'} {'Simple linear'} {'Weighted linear'}]',TAB68);
TAB95 = cat(2,[{'Posterior Kalman (95\%)'} {'Simple linear'} {'Weighted linear'}]',TAB95);
TAB = cat(1,TAB68,TAB95)

TAB_nrmse = [{'NRMSE [\%]'} repmat({' '},1,2)]
TAB_nrmse = cat(1,TAB_nrmse,cat(2,num2cell(nrmse_fit_mean'),repmat({' '},3,2)))
TAB_nrmse = cat(2,[{' '} {'Posterior Kalman'} {'Simple linear'} {'Weighted linear'}]',TAB_nrmse);

tablePoly = l_latexify(TAB,TAB_nrmse,numel(pred_sel),4)
%%
% % prepare latex output
% tablePoly = arr2latex(X*100,{'3.0f'}, ...
%   'collabel',TAB(1,:),...
%   'rowlabel',TAB(2:end,1),...
%   'colspec',colspec, ...
%   'hline','off','centering','on', ...
%   'caption',caption, ...
%   'label',label,...
%   'h2line',3);


savepath = mfilename('fullpath');
savepath = savepath(1:end-19);

tabname = [savepath '../tab/' 'tablePolyfit.tex'];
figname1 = [savepath 'simpleLinear.pdf'];
figname2 = [savepath 'weightedLinear.pdf'];

if save2file
  fileID = fopen(tabname,'w');
  fprintf(fileID,'%s\n',tablePoly);
  fclose(fileID);
  disp(['saved table: ' tabname]);

  figure(1)
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 500 350]);
  print(gcf, figname1, '-dpdf')
  disp(['saved figure: ' figname1]);
  figure(2)
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 500 350]);
  disp(['saved figure: ' figname2]);
  print(gcf, figname2, '-dpdf')
else
  disp(['didn''t save table: ' tabname]);
  disp(['didn''t save figure: ' figname1]);
  disp(['didn''t save figure: ' figname2]);

end



%% ------------------------------------------
function out = l_rolling(vector, window)
% L_ROLLING generates a matrix with the rolling windows
%
% Finds the rolling WINDOW for the data VECTOR. In the output each new row
% is the next roll to consider.
l = length(vector);
m = l - window + 1;
out = vector(hankel(1:m, m:l));

% note however that we could include smaller windows "up to" the start and
% fill the blanks with NaN.
outshorter = NaN(window-1,window);
for k = 1:window-1
  outshorter(k,1:k) = vector(1:k);
end
out = cat(1,outshorter,out);
end

function tex = l_latexify(TAB1,TAB2,N,hmid)
%L_LATEXIFY generates a table environment with a tabular of the
% correct format for the paper result.

    % entries are actually string. We want them as numbers.
    X = cell2mat(TAB1(2:end,2:end));

    %caption = ['RMSE of weekly predictions for proposed model and baseline models'];
    caption = ['Frequency (\%) of all weekly reported predictions for Uppsala' ...
               ' that fell inside of (68/95\% CrI, 7 days ahead) and NRMSE, evaluated on' ...
               ' the following week (N~=~' num2str(N) '), using the' ...
               ' compared models.']

    label = ['tab:modelcomparison'];
    colspec = ['l' 'r' 'r' 'r'];

    tex ='';
    tex = [tex '\begin{table}[htp]' newline ...
           '\centering' newline ...
           '\caption{' caption newline ...
           '2020 to May 2021.}' newline ...
           '\label{' label '}' newline ...
           '\begin{tabular}{' colspec '}'];

    % per row
    TABlist = {TAB1 TAB2};
    for it = 1:numel(TABlist)
        TAB = TABlist{it};
        if it == 1
            tex = [tex '\hline'];
        else
            tex = [tex '\hdashline \hline'];
        end
        for row = 1:size(TAB,1)
            ROW = '';
            % per column, append a '&' sign
            for col = 1:size(TAB,2)
                val = TAB{row,col};
                if isnumeric(val)
                    val = num2str(100*round(val,2,'significant'));
                end



                if col == size(TAB,2) % final column, add '\\' instead of '&'
                    ROW = [ROW val ' \\'];
                    if row == hmid & it == 1
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
    end

    tex = [tex newline ...
      '\hline' newline ...
      '\end{tabular}' newline ...
      '\end{table}'];
end