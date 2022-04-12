%Rposterior creates figure with reproduction number.

% S. Engblom 2021-08-21 (minor revision)
% H. Runvik 2021-08-20

if ~exist('savetofile','var')
    savetofile = true; % print to file false/true
end

% posterior dates to plot
dateend_ = '210531';

shiftstr = '_1';
datestart = 200401; % where to start the plot

regions_ = regions();
% pick one region (if not these specified)
% Uppsala:2, largest:1,10,12; smallest: 8,9,19
for reg = [2]% 1 10 12 8 9 19]
  region = regions_{reg}

  % construct posterior files
  postname = cell(4,1);
  for k = [1 2 3]
    postname{k} = ['slam' dateend_ '_' region '_monthly' ...
      shiftstr '_URDME' num2str(k) '_100.mat'];
    postname{k} = strcat([postpath 'SLAM/perRegion/'],postname{k});
  end

  postname{4} = ['slam' dateend_ '_' region '_monthly' ...
    shiftstr '_100.mat'];
  postname{4} = strcat([postpath 'SLAM/perRegion/'],postname{4});

  R_post_ = [];
  for k = [1 2 3]
    dynoptname_ = ['dynOptPosterior' dateend_ '_' region '_URDME' num2str(k) '.mat'];
    dynoptname_ = strcat([postpath 'dynOpt/'],dynoptname_);
    load(dynoptname_,'R_post')
    R_post_ = cat(2,R_post_,R_post);
  end



  % dynamic beta data
  dynoptname = ['dynOptPosterior' dateend_ '_' region '.mat'];
  dynoptname = strcat([postpath 'dynOpt/'],dynoptname);
  load(dynoptname,'R_post','dates');

  % create a common time frame
  date0 = 200301; % start of common time
  dateend = dates(end); % end of common time from dynOpt
  DATE0 = datenum(2e3+floor(date0/1e4), ...
    mod(floor(date0/1e2),1e2) ,...
    mod(date0,1e2));
  DATEEND = datenum(2e3+floor(dateend/1e4), ...
    mod(floor(dateend/1e2),1e2) ,...
    mod(dateend,1e2));
  DATES = DATE0:DATEEND;
  DATES = datevec(DATES);
  dates = DATES(:,1:3)*[1e4 1e2 1]'-2e7;
  t = 1:numel(dates);

  %
  figure(1), clf, hold on;
  boxw = 15;
  color_softgreen = [64 166 33]/255; % color of line
  color_lightblue = [0.2 0.2 0.8];
  %boxc = color_lightblue; % color of boxes
  boxc = [repmat(color_softgreen,3,1); 0 0 1];
  color_FHMgray = [[120 120 120]/255 0.5];

  % (2) dynamic R
  hline = plot(t(end-numel(R_post)+1:end),R_post,'r');

  %%
%   for k = 1:3
%       hline_ = plot(t(end-numel(R_post_(:,k))+1:end),R_post_(:,k),'Color',boxc(1,:));
%   end
  %%
  R_post_mean = mean(R_post_,2);
  hline_ = plot(t(end-numel(R_post_mean)+1:end),R_post_mean,'--r');

  % (1) the boxplot
  % find centerpoint of each slab
  for k = [1 4]
    if k == 1
      R0 = [];
      for kk = 1:(numel(postname)-1)
        try
          temp = load(postname{kk});
          rates = temp.rates;
          R0 = cat(2,R0,rates.R0);
        end
      end
      rates.R0 = R0; % piggy back on the last rates.
    else
      temp = load(postname{k});
      rates = temp.rates;
    end

    try % temporary fix between legacy files structure.
      slabstop = rates.meta.slabstop;
    catch
      try
        slabstop = temp.slabstop;
      catch
        slabstop = temp.amparam.slabstop;
      end
    end
    if size(slabstop,1) > 1
      slabstop = slabstop(:,1)';
    end

    midpoint = slabstop(1:end-1) + round(diff(slabstop/2));
    if k < 4
      midpoint = midpoint + 1; % rasterize
    end

    if temp.amparam.date(1) == 200319
      xdates = midpoint + 18; % adjusted to include Mars 1 (from 200319)
    elseif temp.amparam.date(1) == 200401
      xdates = midpoint + 31; % as above but from 200401
    else
      error('not defined starting date');
    end



    boxplot(rates.R0','positions',xdates,...
      'widths',boxw,'colors',boxc(k,:),'symbol','','Whisker',1.5);
    if k == 1
      hbox_urdme = findobj(gca,'tag','Box');
    else
      hbox = findobj(gca,'tag','Box');
    end
  end

  %%
% %   % (3) FHM's R-estimates
% %   FHM_R = loadData('FHM_R');
% %
% %   % of FHM_R and C19/RU are of different length. The following code will
% %   % adress that and make sure that the correct data is plotted.
% %   dates_FHM = [find(dates == FHM_R.date(1)) 0];
% %   if dates(end) < FHM_R.date(end)
% %     dates_FHM(2) = numel(dates);
% %     FHM_final = find(FHM_R.date == dates(end));
% %   else
% %     dates_FHM(2) = find(dates == FHM_R.date(end));
% %     FHM_final = numel(FHM_R.date);
% %   end
%
%
%   %dates_FHM = [1 303];
%   hFHM = plot(t(dates_FHM(1):dates_FHM(2)),...
%     FHM_R.mid(1:FHM_final,reg),'color',color_FHMgray);
%   hold off

  xlabel('','Interpreter','Latex');
  ylabel('$R_{t}$','Interpreter','Latex');
  switch reg
    case 10
      region_latex = 'Sk\aa{}ne';
    case 12
      region_latex = 'V\"astra G\"otaland';
    case 19
      region_latex = 'J\"amtland';
    otherwise
      region_latex = region;
  end
  h = title(region_latex,'Interpreter', 'Latex');

  ax = gca;
  ax.TickLabelInterpreter = 'latex';
  yline(1,'--k');

  % for abbreviations, see https://www.bydewey.com/monthdayabb.html
  strmonths = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'June'...
    'July' 'Aug' 'Sept' 'Oct' 'Nov' 'Dec'};
  % 15th each month:
  ixdates = find(mod(dates,100) == 1);
  xxdates = t(ixdates);
  xticks(xxdates);
  % what month:
  months = mod(floor(dates/1e2),1e2);
  slabtitle = strmonths(months(xticks));
  slabtitle(1:2:end) = {''}; % only every 2nd

  % 2021 delimiter
  xline(t(dates == 210101),'--k','2021', ...
    'LabelVerticalAlignment','bottom','LabelOrientation','horizontal', ...
    'interpreter','latex','HandleVisibility','off');

  % xticks & -labels
  xticklabels(slabtitle);
  xtickangle(45);

  % yticks & -labels
  yticklabel = split(num2str(0:0.1:2.4));
  yticklabel(2:2:end) = {''};
  set(gca,'ytick',0:0.1:2.4,'yticklabel',yticklabel);

  % axis fit & legend
  grid on
  if reg == 2 % uppsala can be tighter
    ylim_1 = 0.5;
  else
    ylim_1 = 0.2;
  end
  axis([t(dates == datestart) t(end) ylim_1 2.0])

  %%
  legend([hbox(1) hline(1) hbox_urdme(1) hline_(1)],'Posterior','Quasi-ML', ...
    'Posterior (SB)', 'Quasi-ML (SB)', ...
    'interpreter','latex');
  %%
  % polish target output size
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 500 350]);

  % print to file
  if savetofile
    figure(1)
    printpath = mfilename('fullpath');
    printpath = [printpath(1:end-15) '..' printpath(end-15:end)];
    printpath = [printpath '_' strrep(region,' ','_')];
    printpath = strrep(printpath,'å','a');
    printpath = strrep(printpath,'ä','a');
    printpath = strrep(printpath,'ö','o');
    %print('-depsc',printpath);
    print('-dpdf',printpath);

    disp(['saved figure: ' printpath]);
    % *** does not always work, run manually in this case:
    %unix(['epstopdf ' printpath '.eps']);
  end
end