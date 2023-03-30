% Rposterior creates figure with reproduction number. The figure
% consists of two parts: the posterior per slab in boxplots and the
% marginal boostrap per day as a line. The script requires the
% posterior file(s) to be available otherwise, you user have to see to
% it to generate the needed data.

% S. Engblom 2021-08-21 (minor revision)
% H. Runvik 2021-08-20

if ~exist('savetofile','var')
    savetofile = true; % print to file false/true
end
% pick one region (if not these specified)
% Uppsala:2, largest:1,10,12; smallest: 8,9,19
if ~exist('reg','var')
    reg = [2 1 10 12 8 9 19];
else
    if max(reg) > 21
        error(['Only supported to run 21 regions,'...
               'National posterior is sampled from a basket of region']);
    end
end

if ~exist('verb','var')
    verb=false;
end

if ~exist('dateend_','var')
    dateend_ = 210531;
end
dateend_dyn = dateend_;

if ~exist('interp','var')
    interp=2;
end

if interp == 1
    fileend = '.mat';
elseif interp == 2
    fileend='_update.mat';
else
    error('Only supporting interp 1 or 2')
end
printNtest = true; % write out the number of tests in period 200805-210526

% posterior dates to plot
shiftstr = '_1';
datestart = 200401; % where to start the plot

regions_ = regions(false);
plotid = 0;
for rid = reg
    plotid=plotid+1;
    region = regions_{rid};
    if verb
        disp(['Running: ' region]);
    end

  % construct posterior files
  postname = ['/slam' num2str(dateend_) '_' region '_monthly' ...
              shiftstr '_100' fileend];



  % dynamic beta data
  dynoptname = ['dynOptPosterior' num2str(dateend_dyn) '_' region fileend];
  postname = strcat([postpath 'KLAM/perRegion/'],postname);
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

  % (1) the boxplot
  figure(plotid), clf, hold on;
  boxw = 15;
  color_softgreen = [64 166 33]/255; % color of line
  color_lightblue = [0.2 0.2 0.8];
  %boxc = color_lightblue; % color of boxes
  boxc = [0 0 1];
  color_FHMgray = [[120 120 120]/255 0.5];

  % find centerpoint of each slab
  temp = load(postname);
  rates = temp.rates;
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
  if temp.amparam.date(1) == 200319
    xdates = midpoint + 18; % adjusted to include Mars 1 (from 200319)
  elseif temp.amparam.date(1) == 200401
    xdates = midpoint + 31; % as above but from 200401
  else
    error('not defined starting date');
  end

  boxplot(rates.R0','positions',xdates,...
    'widths',boxw,'colors',boxc,'symbol','','Whisker',1.5);
  hbox = findobj(gca,'tag','Box');

  % (2) dynamic R
  hline = plot(t(end-numel(R_post)+1:end),R_post,'r');

  % (3) FHM's R-estimates
  FHM_R = loadData('FHM_R');

  % of FHM_R and C19/RU are of different length. The following code will
  % adress that and make sure that the correct data is plotted.
  dates_FHM = [find(dates == FHM_R.date(1)) 0];
  if dates(end) < FHM_R.date(end)
    dates_FHM(2) = numel(dates);
    FHM_final = find(FHM_R.date == dates(end));
  else
    dates_FHM(2) = find(dates == FHM_R.date(end));
    FHM_final = numel(FHM_R.date);
  end


  %dates_FHM = [1 303];
  hFHM = plot(t(dates_FHM(1):dates_FHM(2)),...
    FHM_R.mid(1:FHM_final,rid),'color',color_FHMgray);
  hold off

  xlabel('','Interpreter','Latex');
  ylabel('$R_{t}$','Interpreter','Latex');
  switch rid
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
  if interp ==  1
      ylim_2 = 2;
  elseif interp == 2
      if rid == 2
          ylim_2 = 2.4;
      elseif rid == 8
          ylim_2 = 4;
      else
          ylim_2 = 2.5;
      end
  end

  yticklabel = split(num2str(0:0.1:ylim_2));
  yticklabel(2:2:end) = {''};
  set(gca,'ytick',0:0.1:ylim_2,'yticklabel',yticklabel);

  % axis fit & legend
  grid on
  ylim_1 = 0.2;



  axis([t(dates == datestart) t(end) ylim_1 ylim_2])
  legend([hbox(1) hline hFHM],'Posterior','Quasi-ML', ...
    'PHA','interpreter','latex');

  % polish target output size
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 500 350]);

  % print to file
  if savetofile
    figure(plotid)
    printpath = mfilename('fullpath');
    printpath = [printpath(1:end-10) '..' printpath(end-10:end)];
    printpath = [printpath '_' strrep(region,' ','_')];
    printpath = strrep(printpath,'å','a');
    printpath = strrep(printpath,'ä','a');
    printpath = strrep(printpath,'ö','o');
    print('-dpdf',printpath);

    if verb
        disp(['saved figure: ' printpath])
    end
  end

end
%%
if printNtest
  load 'Ncounties'
  N100k = sum(N,1)/1e5;
  test = loadData('FHM_test');
  tspan_test = find(test.date == 200805):find(test.date == 210526);
  Ntest = sum(test.Nper100k(tspan_test,reg).*N100k(reg),1);
  Ntest_c = sprintfc('%d',round(Ntest,4,'significant')');
  Ntest100k = sum(test.Nper100k(tspan_test,reg),1);
  Ntest100k_c = sprintfc('%d',round(Ntest100k,5,'significant')');
  names = regions();
  Ntest_tab = cat(2,names(reg),Ntest_c, Ntest100k_c);
  if verb
      disp(Ntest_tab);
  end
end
