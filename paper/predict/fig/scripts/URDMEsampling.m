%URDMESAMPLING URDME synthetic data for all 21 regions at once.


% R. Marin 2022-03-25 (adjusting for plotting)
% S. Engblom 2021-05-28

% change these for saving, replicas, and regions
if ~exist('savetofile','var')
    savetofile = false;
end

includeUDS = false; % don't generate samples using the uds solver
if ~exist('Nreplicas','var')
    Nreplicas = 25;
end


if ~exist('reg','var')
    reg = 1:21;
end

if ~exist('regplot','var')
    regplot = [2 1 10 12 8 9 19];
end

if ~all(ismember(regplot,reg))
    error('all regions that you intend to plot needs to be generated as well');
end


% note: comment away 'solver', 'reg', and consider to pre-compile once
% (set 'compile' = 0 when conctructing the umod struct)
solver = 'ssa';
date = '210531';
stopdate = '210531';
divide = true; % I -> H, W -> H2, D = {D_I, D_H, D_W};
Psamples = inf;% % inf: mean posterior, else: N # of posterior samples.

% construct the rates to run
P = load('slam210531_mean_monthly_1');
P = P.rates;
DR = load('dynOptPosterior210531_all');

% check: that we are talking about the same thing:
if any(P.meta.hash ~= DR.postrates.meta.hash) || ...
      any(P.meta.dataHash ~= DR.postrates.meta.dataHash)
  error('Posterior mismatch.');
end

% get daily beta
denum = (DR.postrates.thetaE_/DR.postrates.sigma+ ...
         (DR.postrates.F0ave+(1-DR.postrates.F0ave).*DR.postrates.F1ave)./DR.postrates.gammaI+ ...
         (1-DR.postrates.F0ave).*DR.postrates.thetaA_./DR.postrates.gammaA);
R = DR.R_all';
beta = R.^2./denum(DR.postrates.meta.slabs);
% *** can be solved by getC19beta/dailyrates (untested & unreviewed)

u0 = DR.xSim_all(reg,:,1)';
phi0 = DR.xSim_all(reg,4,1)';
if divide
  u0 = [u0; zeros(3,numel(reg))];
end


% main call
disp(' generating data ...');
disp('... ssa samples');
covid19enger_run_post

%
% pack up output
clear D
D.regions = regions(reg);
D.date = DATES;
D.U = reshape(permute(umod.U(:,:,:,:),[2 1 3 4]),numel(DATES),Nspecies,Nvoxels,Nreplicas,[]);
D.vars = umod.private.Species;

% save rates
Rates.beta = beta;
Rates.R = R;
if Psamples > 1
  Rates.sigma = P.sigma;
  Rates.gammaI = P.gammaI;
  Rates.gammaH = P.gammaH;
  Rates.gammaW = P.gammaW;

  Rates.F0 = P.F0ave;
  Rates.F1 = P.F1ave;
  Rates.F2 = P.F2ave;
  Rates.F2d = P.F2dave;
  Rates.F3 = P.F3ave;
  Rates.F3d = P.F3dave;
  Rates.F4 = P.F4ave;
end
Rates.rho = P.rho;
Rates.thetaI = P.thetaI;
Rates.thetaA = P.thetaA;
Rates.thetaE = P.thetaE;
Rates.thetaI_ = P.thetaI_;
Rates.thetaA_ = P.thetaA_;
Rates.thetaE_ = P.thetaE_;
Rates.lambda = P.lambda;
D.rates = Rates;

% .meta-info
D.meta = P.meta;
D.meta.posteriorHash = P.meta.hash;
D.meta.hash = fsetop('check',D.U(:));

save('URDME_all','D');
%
if includeUDS
    disp('... mean field approximation');
    % run Euler forward
    solver = 'uds';
    temp = Nreplicas; % we only need 1 sample from UDS.
    Nreplicas = 1;
    covid19enger_run_post
    Nreplicas = temp;
    load('URDME_all','D');
    D.EulFwd = reshape(permute(umod.U,[2 1 3]),numel(DATES),Nspecies,Nvoxels,Nreplicas);
    Nspecies_ = 7;
    D.dynOpt = reshape(permute(DR.xSim_all(reg,:,:),[3 2 1]),numel(DATES),Nspecies_,Nvoxels,Nreplicas);

    save('URDME_all','D');
end

%%
% plot simulations with data
disp('generating plots ...');

if divide
  Hcol = [5 9];
else
  Hcol = [5];
end
Wcol = [6];

% if numel(reg) == 21
%   tiledlayout(7,3)
% else
%   tiledlayout(3,3)
% end

alpha=0.1;
for regid=regplot
  %nexttile;
  H = squeeze(sum(D.U(:,Hcol,regid,:),[2 3]));
  W = squeeze(sum(D.U(:,Wcol,regid,:),[2 3]));
  plot(H,'Color',[0 0 1 alpha],'HandleVisibility','off'), hold on
  plot(W,'Color',[1 0 0 alpha],'HandleVisibility','off')

  if includeUDS
      Hmean = sum(D.EulFwd(:,Hcol,regid),[2 3]);
      Wmean = sum(D.EulFwd(:,Wcol,regid),[2 3]);
  else
      Hmean = mean(sum(D.U(:,Hcol,regid,:),[2 3]),4);
      Wmean = mean(sum(D.U(:,Wcol,regid,:),[2 3]),4);
  end
  plot(Hmean,'b')
  plot(Wmean,'r')

  data = loadData('C19');
  tspan_post = find(data.date == D.date(1)):find(data.date == D.date(end));
  plot(sum(data.H(tspan_post,reg(regid)),2),'.b')
  plot(sum(data.W(tspan_post,reg(regid)),2),'.r')
  title(l_getlatex(reg(regid)),'Interpreter', 'Latex')

  % Apply correct formating
  ax = gca;
  ax.TickLabelInterpreter = 'latex';

  x0 = 1:numel(tspan_post);
  xtk = fliplr(x0(end:-28:1));
  %xtk = fliplr(tspan_post(end:-28:1));

  set(gca,'xtick',xtk);


  % wanted xticks
  % for abbreviations, see https://www.bydewey.com/monthdayabb.html
  strmonths = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'June'...
    'July' 'Aug' 'Sept' 'Oct' 'Nov' 'Dec'};

  % 15th each month:

  ixdates = find(mod(D.date,100) == 1);
  tspan_urdme = 1:numel(D.date);
  xxdates = tspan_urdme(ixdates);
  dates=D.date(xxdates);
  xticks(xxdates);
  % what month:
  months = mod(floor(dates/1e2),1e2);
  slabtitle = strmonths(months);
  slabtitle(2:2:end) = {''}; % only every 2nd

  % xticks & -labels
  xticklabels(slabtitle);
  xtickangle(45);


  % 2021 delimiter
  xline(find(D.date == 210101),'--k','2021', ...
    'LabelVerticalAlignment','top','LabelOrientation','horizontal', ...
    'interpreter','latex','HandleVisibility','off');


  legentries = {'Hospital [H]','Intensive [W]'};
  leg= legend(legentries, ...
    'Location','northwest','interpreter','latex','FontSize',12,...
    'Orientation','horizontal',...
    'NumColumns', 1);



  grid on
  axis tight

  switch regions{regid}
    case 'Gotland'
      ylim([0,25]);
    case 'Blekinge'
      ylim([0,60]);
    case 'Jämtland'
      ylim([0,35]);
  end

  hold off


  %
  % ***************************************
  % *** save to files ***
  % ***************************************
  savepath = mfilename('fullpath');
  savepath = [savepath(1:end-21) 'URDME_samples'];
  savepath = [savepath '_' strrep(regions{reg(regid)},' ','_')];
  savepath = strrep(savepath,'å','a');
  savepath = strrep(savepath,'ä','a');
  savepath = strrep(savepath,'ö','o');


  if savetofile % true if we should compute the error table




    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Position',[100 100 500 350]);


    print('-dpdf', savepath)
    disp(['saved figure: ' savepath]);


  else
    disp(['did not save figure:' savepath])
  end

end
disp(' ... done with plots');







function latexname = l_getlatex(regid)
%L_GETLATEX quick and dirty to get the translation from text to latex
%  friendly Swedish's region names.
regionList_tex = {'Stockholm' 'Uppsala' 'S\"{o}dermanland' '\"{O}sterg\"{o}tland' ...
  'J\"{o}nk\"{o}ping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
  'Sk\aa{}ne' 'Halland' 'V\"{a}stra G\"{o}taland' 'V\"{a}rmland' '\"{O}rebro' ...
  'V\"{a}stmanland' 'Dalarna' 'G\"{a}vleborg' 'V\"{a}sternorrland' ...
  'J\"{a}mtland' 'V\"{a}sterbotten' 'Norrbotten'}';
latexname = regionList_tex{regid};
end