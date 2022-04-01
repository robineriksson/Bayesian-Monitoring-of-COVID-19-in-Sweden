%WEEKLY_PREDICTION Script for generating predictions.

Nprior = 1e3;
plotfun = @plot;
savetofile = true;

 abspath = mfilename('fullpath');
if ~exist([abspath(1:end-9) 'priorpred.mat'], 'file')
  disp('generating priorpred.mat')
  

  lag = 7;

  
  % use CSSS?
  useCSSS=false;
  
  
  % for reference
  region = 'Uppsala'
  
  % posterior file date
  posteriordate = '210531';
  ending         = '1_100'; % at what date does the slabs start.
  
  % dates that "date" the data used
  datadate     = [200401 210531]; % change end date here
  
  
  datasource = 'RU';
  
  
  % folder for posteriors
 
  prefix = [postpath 'SLAM/'];
  
  % Population data
  load Ncounties
  Npop = sum(N,1);
  
  
  % convert region into "test" for data extraction, later.
  regionList = regions();
  test = find(strcmp(region,regionList)) + 3;
  if isempty(test)
    if strcmp(region,'Sweden')
      test = 0;
      
      % construct the Sweden posteriorfile
      if ~exist([prefix 'slam' posteriordate '_Sweden_monthly_' ending '.mat'],'file')
        regs = 1:21;
        files = cell(size(regs));
        for r = 1:numel(regs)
          files{r} = [prefix 'perRegion/slam' posteriordate '_' regionList{regs(r)} '_monthly_' ending '.mat'];
        end
        
        Weights = slweight(files);
        rates = posteriorenger(100,files,Weights);
        save([prefix 'slam' posteriordate '_Sweden_monthly_' ending '.mat'],'rates')
        
        
      end
      
      
    else
      error(['Missing region: ' region]);
    end
  end
  
  
  if test > 0
    posterior = 'perRegion/';
  else
    posterior = '';
  end
  
  %posterior = [posterior '/smc' posteriordate '_22_' region '_monthly'];
  %posterior = [posterior '/slam' posteriordate '_22_' region '_monthly'];
  posterior = [posterior 'slam' posteriordate '_' region '_monthly'];
  
  if useCSSS
    posterior = [posterior '_csss'];
  end
  
  posterior = [posterior '_' ending '.mat'];
  
  % sample rates
  rng(0);
  rates = posteriorenger([],[prefix posterior]);
  
  try
    slabs = rates.meta.slab';
  catch
    slabs = rates.meta.slabs;
  end
  if strcmp(region,'Sweden')
    slabs = slabs(1,:); % assume all the regions, had the same slabs.
  end
  
  
  rates = priorenger(Nprior,true,numel(rates.meta.slabstop)-1);
  rates.meta.date = datadate;
  posterior = 'perRegion/prior_Uppsala.mat';
  % transmission matrix D
  load Dcounties
  
  
  test_ = test; % for plotting keep test_
  
  % load filter data
  if strcmp(datasource,'URDME') % synthetic data
    filename = mfilename('fullpath');
    Data_raw = load([filename(1:end-24) 'data/sources/URDME/' 'URDME210418_Stockholm_Uppsala.mat']);
    
    Data = struct();
    Data.regions = Data_raw.D.region;
    Data.date = Data_raw.D.date;
    Data.H = Data_raw.D.U(:,5,1);
    Data.W = Data_raw.D.U(:,6,1);
    Data.D = Data_raw.D.U(:,7,1);
    
    Data.I = Data_raw.D.U(:,1,1); % currently for plotting.
    
    Data.hash = -1;
    Data.rev = -1;
    Data.reg = 'URDME210418_Stockholm_Uppsala';
    T = 1;
    D = 1;
    lan = {Data.regions};
    test = 0;
  else
    Data = loadData(datasource);
    Data = polishData(Data,'D','Dinc',1);
    Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});
    
    % this is to allow for perRegion posteriors
    if test >= 4
      D = 1;
      lan = 1;
      Data.D = Data.D(:,strcmp(Data.regions, region));
      Data.W = Data.W(:,strcmp(Data.regions, region));
      Data.H = Data.H(:,strcmp(Data.regions, region));
      Npop = Npop(strcmp(Data.regions, region));
      Data.regions = Data.regions(test-3);
      test = 0;
    end
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
  
  % *** Code to allow for CSSS data, a bit measy, should be looked into.
  if useCSSS
    % correction factor in xmod = fac*x
    load Fcases_interp
    fac = sum(Eave(3:end))/sum(Iave(3:end));
    
    % CSSS-data
    if strcmp(region,'Uppsala')
      csss = loadData('CSSS_RU');
      Ihi = csss.Ihi(:,1);
      Ilo = csss.Ilo(:,1);
      Imid = csss.Imid(:,1);
    else
      csss = loadData('CSSS');
      Ihi = csss.Ihi(:,strcmp(Data.regions, region));
      Ilo = csss.Ilo(:,strcmp(Data.regions, region));
      Imid = csss.Imid(:,strcmp(Data.regions, region));
    end
    ixcsss = find(DATES == csss.date(1)):find(DATES == csss.date(end));
    tspan_csss = TSPAN(ixcsss);
    
    %   Ihi = csss.Ihi;
    %   Ilo = csss.Ilo;
    %   Imid = csss.Imid;
    
    Imid(csss.date > 201231) = NaN; % data post 31 Dec 2020, is "bad".
    
    % settle on a relative variance from given asymmetric 95% CI
    varI = mean(((Ihi-Ilo)/2./Imid).^2); % generally larger
    %varI = max(varI);
    %varI = varI(test-3);
    %   varI(isnan(varI)) = [];
    %   varI = mean(varI);
    %   %varI = mean(((Ihi-Ilo)/4./Imid))^2;
    
    % specify observation model
    obsrates.states = {5 6 7}; % state #5, #6, #7
    obsrates.indobs = {1}; % CSSS data for the I-compartment
    obsrates.indobspars = {100./(sum(Npop)*fac)}; % probably WRONG?!
    obsrates.nstate = 8;
    obsrates.R0 = [ones(3,1); 1e-4];
    obsrates.rdiag = [(1e-3)^2*ones(3,1); varI];
    
    % insert CSSS data
    Idata = permute([NaN(ixcsss(1)-ixdata(1),numel(lan)); Imid; ...
      NaN(ixdata(end)-ixcsss(end),numel(lan))],[3 2 1]);
    Idata_hi = permute([NaN(ixcsss(1)-ixdata(1),numel(lan)); Ihi; ...
      NaN(ixdata(end)-ixcsss(end),numel(lan))],[3 2 1]);
    Idata_lo = permute([NaN(ixcsss(1)-ixdata(1),numel(lan)); Ilo; ...
      NaN(ixdata(end)-ixcsss(end),numel(lan))],[3 2 1]);
    Ydata = cat(1,Ydata,Idata);
  else % don't use CSSS.
    % specify observation model
    obsrates.states = {5 6 7}; % state #5, #6, #7
    obsrates.indobs = {}; % No indirect measurements
    obsrates.indobspars = {};
    obsrates.nstate = 8;
    obsrates.R0 = 1*ones(3,1);
    obsrates.rdiag = (1e-3)^2*ones(3,1);
  end
  % specify output model
  H = getC19obs(obsrates); % i.e., "H", "W", "D" ("I" if CSSS)
  % read as: "Stockholm", "Uppsala", (...), "Sweden total"
  %T = sparse([1 2 3*ones(1,numel(lan))],[1 2 1:numel(lan)],ones(1,2+numel(lan)));
  T = [speye(numel(lan)); sparse(1,1:numel(lan),1)];
  G = kron(T,H);
  
  % uncertainty parameters
  Q = struct();
  Q.Q0 = speye(numel(lan)*obsrates.nstate);
  Q.qdiag = 0.05^2;
  
  
  % Added rows for outputting phi
  %if phiplotting
  % order follows {I A E phi H W D R}
  Gphi = kron(speye(numel(lan)),[0 0 0 1 0 0 0 0]);
  G = [G; Gphi];
  %end
  
  % Added rows for outputting R
  %if Rplotting
  GR = kron(speye(numel(lan)),[0 0 0 0 0 0 0 1]);
  G = [G; GR];
  %end
  
  %if Iplotting && ~useCSSS
  GI = kron(speye(numel(lan)),[1 0 0 0 0 0 0 0]);
  G = [G; GI];
  %end
  
  % these were used durring inference (some regions are sensitive)
  exception.LB = -1e2;
  exception.UB = 1e7;
  exception.SDFAC = 0.25;
  exception.AbsMagn = 1e4;
  
  % with historic lag or not
  
  disp(['qdiag: ' num2str(Q.qdiag)]);
  [Z_,covZ_,Z,covZ] = ...
    C19filt_lag(G,rates,D,obsrates,Ydata,slabs,...
    numel(ixfilter),lag,Q);
  covZ.covZ = []; % save some space
  covZ_.covZ = [];
  
  % when loading, are we loading the posterior we think we're doing?
  meta = struct();
  meta.postHash = rates.meta.hash; % same posterior
  meta.dataHash = Data.hash;  % same data
  
  save([abspath(1:end-9) 'priorpred.mat']);
  disp('saved priorpred.mat')
  
else
  load([abspath(1:end-9) 'priorpred.mat'])
  disp('loading priorpred.mat')
end

%if Iplotting
% Separating I from the rest
ZI = Z(end-numel(lan)+1:end,:,:);
stdZI = covZ.stdZ(end-numel(lan)+1:end,:,:);
Z(end-numel(lan)+1:end,:,:) = [];
covZ.stdZ(end-numel(lan)+1:end,:,:) = [];
Z_(end-numel(lan)+1:end,:,:) = [];
covZ_.stdZ(end-numel(lan)+1:end,:,:) = [];
%end

%if Rplotting
% separating R from the rest
ZR = Z(end-numel(lan)+1:end,:,:);
stdZR = covZ.stdZ(end-numel(lan)+1:end,:,:);
Z(end-numel(lan)+1:end,:,:) = [];
covZ.stdZ(end-numel(lan)+1:end,:,:) = [];
Z_(end-numel(lan)+1:end,:,:) = [];
covZ_.stdZ(end-numel(lan)+1:end,:,:) = [];
%end
%if phiplotting
% separating phi from the rest
Zphi = Z(end-numel(lan)+1:end,:,:);
stdZphi = covZ.stdZ(end-numel(lan)+1:end,:,:);
Z(end-numel(lan)+1:end,:,:) = [];
covZ.stdZ(end-numel(lan)+1:end,:,:) = [];
Z_(end-numel(lan)+1:end,:,:) = [];
covZ_.stdZ(end-numel(lan)+1:end,:,:) = [];
%end

%%
figure(1), clf;
stdZ = covZ.stdZ; % for brevity

% ***************************
% *** Postprocessing code ***
% ***************************

Hcomp = 1; Wcomp = 2; Dcomp = 3;

Icomp= [];

% Full simulation (Sweden) or single region (no transport).
Ydata_swe = squeeze(sum(Ydata,2));

% summed counties
end_ = size(Z,1);

ind = end_-2:end_;
Y_swe = Z(ind,:,:);
stdY_swe = stdZ(ind,:,:);

% scaling for plots
Ydata_swe_plot = Ydata_swe;
smoothing = @(x) sgolayfilt(x,0,7);
Ydata_swe_plot(Dcomp,:) = cat(2,ones(1,1), diff(Ydata_swe(Dcomp,:),1,2));
Ydata_swe_plot(Dcomp,:) = smoothing(Ydata_swe_plot(Dcomp,:));

Y_swe_plot = Y_swe;

Y_swe_plot(Dcomp,:,:) = cat(2,ones(1,1,size(Y_swe,3)), diff(Y_swe(Dcomp,:,:),1,2));
Y_swe_plot(Dcomp,:,:) = smoothing(Y_swe_plot(Dcomp,:,:));



stdY_swe_plot = stdY_swe;



if strcmp(datasource,'URDME')
  Ydata_swe_plot(Icomp,:) = Data.I;
end
filtval = 0;
Y_sweF = max(Y_swe,filtval); % filtered values
Y_sweF_plot = max(Y_swe_plot,filtval); % filtered values



% number of samples to draw (illustration only)
Nsamples_ = 100;


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
plotfun(tspan_filter,squeeze(Y_swe_plot(Dcomp,:,1:Nsamples_)),'-','Color',[color(3,:) 0.1], ...
  'HandleVisibility','off')


% sample predictive mean

comp_vec = [Hcomp Wcomp Dcomp Icomp];

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
hold off,

  xlim(tspan_data([1 end]))
  xtk = fliplr(tspan_data(end:-28:1));
  yyaxis right % no negative axis
  yvr = ylim;
  ylim([0 yvr(2)])
  
  yyaxis left
  ylim([0,400])
  
  



% 
% % below is a trick to make the grid align both with the left and right
% % yaxis.
% 
% yyaxis left
% ytl = get(gca, 'YTick'); % Get Controlling Left Ticks
% ytlv = round(linspace(min(ylim),max(ylim),numel(ytl)));
% ytlc = compose('%d',ytlv); % Tick Label Cell Array
% set(gca, 'YTick',ytl, 'YTickLabel',ytlc)  % Set New right Tick Labels
% 
% 
% 
% yyaxis right
% ylim0 = ylim;
% ytr = get(gca, 'YTick'); % Get Right Tick Values
% 
% tinc = round(diff(ylim)/numel(ytl));  % Tick Increment
% 
% 
% yyaxis right
% tinc0 = round(diff(ylim)/numel(ytl));  % Tick Increment
% ytrv = min(ylim):tinc0:(min(ylim) + tinc0*(numel(ytl)-1));
% 
% 
% 
% ylim([min(ytrv) max(ytrv)])% Set New ‘ylim’
% ytrc = compose('%d',ytrv); % Tick Label Cell Array
% set(gca, 'YTick',ytrv, 'YTickLabel',ytrc)  % Set New Right Tick Labels
% grid on
% 
% if min(ylim) + tinc0*(numel(ytl)-1) < max(ylim0)
%   % top is missing, add one tick on left and right.
%   
%   
%   yyaxis left
%   ytl = get(gca, 'YTick'); % Get Controlling Left Ticks
%   tinc = diff(ytl(1:2));
%   ytlv = [ytl ytl(end) + tinc];
%   ytlc = compose('%d',ytlv); % Tick Label Cell Array
%   set(gca, 'YLim', [min(ytlv) max(ytlv)])
%   set(gca, 'YTick',ytlv, 'YTickLabel',ytlc)  % Set New left Tick Labels
%   
%   
%   yyaxis right
%   ytr = get(gca, 'YTick'); % Get Controlling Left Ticks
%   tinc = diff(ytr(1:2));
%   ytrv = [ytr ytr(end) + tinc];
%   ytrc = compose('%d',ytrv); % Tick Label Cell Array
%   set(gca, 'YLim', [min(ytrv) max(ytrv)])
%   set(gca, 'YTick',ytrv, 'YTickLabel',ytrc)  % Set New Right Tick Labels   end
% end
% % end of trick.


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



leg = {'Hospital [H]','Intensive [W]','Dead/day [D]'};

numcols = 3;


legend(leg, ...
  'Location','northwest','interpreter','latex','FontSize',12,...
  'Orientation','horizontal',...
  'NumColumns', 1);




hold off



savepath = abspath(1:end-17);
if savetofile % true if we should compute the error table
  
  
  
  % finalize the plot
  h = gcf;
  ax = gca;
  ax.TickLabelInterpreter = 'latex';
  %      set(h,'PaperOrientation','landscape');
  %      set(h,'PaperUnits','normalized');
  %      set(h,'PaperPosition',[0 0 1 1]);
  %      set(h,'PaperUnits','centimeter')
  %      set(h,'PaperSize',[20 16])
  %
  % polish target output size
  set(h,'PaperPositionMode','auto');
  set(h,'Position',[100 100 500 350]);
  
  %print(h, ['~/Desktop/lag_' region '.eps'], '-depsc')
  print(h, [savepath 'priorpred.pdf'], '-dpdf')
  
  %
end