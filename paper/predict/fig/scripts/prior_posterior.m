% PRIOR_POSTERIOR plots of marginal densities, per single or
% collection of regions. If a collection, say all regions, is given
% (REG = [...]), the posterior output is weighted according to
% population size. Included in the generated figure is also the prior
% and the mean of the prior and posterior.
%

% S. Engblom 2021-04-03 (Revision)
% S. Engblom 2020-10-26 (Minor revision)
% S. Engblom 2020-10-24


if ~exist('reg','var')
    reg = 1:21;
end

if ~exist('savetofile','var')
    savetofile=false;
end


% load posterior
Nsample = 1e4;
path = postpath;
prefix = 'KLAM/perRegion/';
% "super-posterior"
date = '210531';
ending = '_1_100';
regionList = regions(false);

includeBias=false; % true if to include the bias estimate.

regname = regionList(reg);
files = strcat('slam',date,'_',regname,'_monthly', ending, '.mat');
load Ncounties;


filename = strcat(path,prefix,files);

Weights = sum(N(:,reg));
ratesPosterior = posteriorenger(Nsample,filename,Weights,true);
if includeBias
    load([postpath 'SLAM/marginal_bias'],'marginal');
end
% filetrue = [postpath 'SLAM/slam210531_mean_monthly' ending];
% rates=ratesPosterior;
% save(filetrue,'rates')
meansPosterior = struct;

% standard use of a single posterior (but one can use the above
% code more generally with, e.g., reg = 2!)
%postname = [path(1:end-32) prefix 'slam210328_22_Upp
%postname = 'smc210420_22_Stockholm_monthly_15.mat';
%postname = 'slam210425_22_Stockholm_monthly_22.mat';
%postname = 'slam210418_22_Stockholm_monthly_15.mat';
%postname = 'smc210420_22_Uppsala_monthly_15.mat';
%postname = 'slam210418_22_Uppsala_monthly_15.mat';
%postname = 'slam210418_22_Stockholm_monthly_15.mat';
%postname = 'slam210412_22_Uppsala_monthly_8.mat';
%postname = 'slam210411_22_Stockholm_monthly_8.mat';
%postname = 'slam210404_22_Uppsala_monthly.mat';
%postname = 'slam210404_22_Stockholm_monthly.mat';
%postname = 'slam210405_22_Uppsala_monthly.mat';
%postname = 'slam210405_22_Uppsala_monthly_csss_fix.mat';
%postname = 'slam210405_22_Uppsala_monthly_csss.mat';
% $$$ filename = [path(1:end-32) prefix postname];
% $$$ ratesPosterior = posteriorenger(Nsample,filename);
% $$$ meansPosterior = struct; % built on the fly
%%
% load prior
fix = false;
slabs = ratesPosterior.meta.slabs(1,:);
nslab = max(slabs);
rng(0)

try
  hypfile = ratesPosterior.meta.hypfile;
  if contains(hypfile,'+') % if combined posterior
    stop = strfind(ratesPosterior.meta.hypfile,'+')-2;
    stop = stop(1);
    hypfile = extractBetween(ratesPosterior.meta.hypfile,1,stop);
    hypfile = hypfile{:};
  end

catch
  hypfile = [];
end
ratesPrior = priorenger(Nsample,fix,nslab,hypfile);
meansPrior = priorenger(inf,fix,nslab,hypfile);

% selection of parameters to display
rateNamesSubset =  {'sigma' 'gammaI' 'gammaH' 'gammaW' ...
                    'R0' 'thetaA_' 'thetaE_' ...
                    'half_life' 'E2I' 'HOSP' 'IC_HOSP' 'IFR'};
rateNamesSubset_tex = {'$\sigma^{-1}$ [days]' ...
                      '$\gamma_I^{-1}$ [days]' ...
                      '$\gamma_H^{-1}$ [days]'...
                      '$\gamma_W^{-1}$ [days]' ...
                      '$R_t$' ...
                      '$\theta_A$ [$A \rightarrow \varphi$]' ...
                      '$\theta_E$ [$E \rightarrow \varphi$]' ...
                      'Infectious half-life [hrs]' ...
                      '$E_2I$ [$E \rightarrow I$, \%]' ...
                      'HOSP [$I \rightarrow H$, \%]' ...
                      'IC HOSP [$H \rightarrow W$, \%]' ...
                      'IFR [$E \rightarrow D$, \%]'};

figure(1), clf,
for i = 1:numel(rateNamesSubset)
    ha(i) = subplot(3,4,i); hold on,
    set(0,'defaultfigurecolor',[1 1 1]);
    set(0,'defaulttextinterpreter','latex');
    set(0,'DefaultTextFontname','CMU Serif');
    set(0,'DefaultAxesFontName','CMU Serif');
    set(gca,'TickLabelInterpreter','latex');

    name = rateNamesSubset{i};
    prior = ratesPrior.(name)(end,:);
    prmean = meansPrior.(name)(end);
    posterior = mean(ratesPosterior.(name),1);
    pomean = mean(posterior);
    if includeBias
        pobias = sqrt(mean(cell2mat(marginal(fsetop('ismember',marginal(:,2,1),{name}),1,reg))));
        if any(strcmp({'R0','IFR'},name))
            pobias = pobias(end);
        end
    end

    switch name
      case {'F0' 'HOSP' 'IC_HOSP' 'IFR' 'E2I'}
        prior = prior.*100;
        prmean = prmean*100;
        posterior = posterior.*100;
        pomean = pomean*100;
        if includeBias, pobias = pobias*100; end,
      case {'sigma' 'gammaI' 'gammaH' 'gammaW'}
        prior = 1./prior;
        prmean = 1/prmean;
        posterior = 1./posterior;
        pomean = mean(posterior); % (need not commute)
        if includeBias
            % if inverse, we need to 'correct' the bias
            temp = 1/pomean + pobias;
            pobias = abs(1/temp - pomean);
        end

      case 'half_life'
        prior = prior*24;
        prmean = prmean*24;
        posterior = posterior*24;
        pomean = pomean*24;
        if includeBias
            pobias= pobias*24;
        end
    end

%     disp(name)
%     disp(quantile(posterior,[0.025,0.975]))
%     disp(pomean)

    % posterior
    [fp,xip] = ksdensity(posterior,'Function','pdf');
    xx = [xip,flip(xip,1)];
    ff = [fp,flip(fp,1)];
    h_post = fill(xx,ff,[0 0 1]);
    hold on
    set(h_post,'facealpha',.1)

    % prior
    [f,xi] = ksdensity(prior,'Function','pdf');
    xx = [xi flip(xi,1)];
    ff = [f flip(f,1)];
    h_prior = fill(xx,ff,[1 0 0]);
    set(h_prior,'facealpha',.1)

    xlabel(rateNamesSubset_tex{i},'interpreter','latex');

    % mean vertical line, with annotation.
    hline_prior = xline(prmean,'--','LineWidth',1,'Color', [1 0.2 0.2]);
    hline_post = xline(pomean,'-','LineWidth',2,'Color', [0 0.4470 0.7410]);
    str = ['$=$ ' num2str(round(pomean,2,'significant'))];
    if includeBias
        str = [str '$\pm$' num2str(round(pobias,2,'significant'))];
    end
    b(i) = annotation('textbox','String',str,'Position',ha(i).Position,...
      'Vert','top','HorizontalAlignment','right',...
      'FitBoxToText','on','interpreter','latex',...
      'BackgroundColor',[1 1 1], ...
      'FaceAlpha',0.01,...
      'EdgeColor','None','Margin',0,'FontSize',8,...
      'Color',[0 0.4470 0.7410]);
    set(gcf,'color','w');

    % minor improvements on axis
    switch name
      case 'R0'
        xlim([0.9 1.1]);
        str = {['E[prior]=' num2str(round(prmean,2,'significant'))] ''};

        b(i) = annotation('textbox','String',str,'Position',ha(i).Position,...
          'Vert','bot','HorizontalAlignment','right',...
          'FitBoxToText','on','interpreter','latex',...
          'BackgroundColor',[1 1 1], ...
          'FaceAlpha',0.01,...
          'EdgeColor', 'None','Margin',2,'FontSize',7,...
          'Color', [1 0.2 0.2]);
      case 'HOSP', xlim([0 10])
      case 'IC_HOSP', xlim([0 40])
    end
    set(gca,'ytick',[]);
    set(get(gca,'Xruler'),'Exponent',0)
    hold off


end

legend([h_post(1), h_prior(1)], ...
  'Posterior','Prior',...
  'Orientation','horizontal','interpreter','latex', ...
  'Position',[0.47 0.95 0.1 0.03],'FontSize',14);
%
% legend('Posterior','Prior', ...
%   'Orientation','horizontal','interpreter','latex', ...
%   'Position',[0.47 0.95 0.1 0.03],'FontSize',14);
legend boxoff

% polish target output size
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 700 400]);

printpath = mfilename('fullpath');
printpath = [printpath(1:end-23) 'posterior_sweden.pdf'];
% print to file
if savetofile
  %print('-depsc',printpath);
  print('-dpdf',printpath);
  % *** does not always work, run manually in this case:
  %unix(['epstopdf ' printpath '.eps']);
  disp(['saved figure: ' printpath]);
else
  disp(['didn''t save figure: ' printpath]);
end