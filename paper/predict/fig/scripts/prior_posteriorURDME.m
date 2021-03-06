%PRIOR_POSTERIORURDME generate figure that compares the actual data posterior with
% the ones recovered from the URDME equivalents. The script is a derivative
% of PRIOR_POSTERIOR as it only adds the URDME posterior on top of the
% existingposterior plot.

% R. Eriksson 2022-01-27

if ~exist('reg','var')
    reg = 1:21;
else
    if max(reg) > 21
        error(['Only supported to run 21 regions,'...
               'National posterior is sampled from a basket of region']);
    end
end

if ~exist('savetofile','var')
    savetofile=false;
end

if ~exist('verb','var')
    verb=false;
end

% load posterior
Nsample = 1e4;
path = postpath;
prefix = 'KLAM/perRegion/';
% "super-posterior"
date = '210531';
ending = '_1';
ending_ = '_100';
ending = [ending ending_];
regionList = regions(false);


regname = regionList(reg);
files = strcat('slam',date,'_',regname,'_monthly', ending, '.mat');

% pull all bootstrapped samples
files_urdme = cell(0,1);
for k = 1:3

    ending_urdme = ['_1_URDME',num2str(k) ending_];
    if strcmp(ending_,''), prefix_urdme = 'URDME/'; else, prefix_urdme = ''; end
    files_urdme = cat(1,files_urdme,strcat(prefix_urdme,'slam',date,'_',regname,'_monthly', ending_urdme, '.mat'));
end

load Ncounties;


filename = strcat(path,prefix,files);
filename_urdme = strcat(path,prefix,files_urdme);
Weights = sum(N(:,reg));
ratesPosterior = posteriorenger(Nsample,filename,Weights,true);
ratesPosterior_urdme = posteriorenger(Nsample,filename_urdme,repmat(Weights,1,3),true);
% filetrue = [postpath 'SLAM/slam210531_mean_monthly' ending];
% rates=ratesPosterior;
% save(filetrue,'rates')
meansPosterior = struct;
meansPosterior_urdme = struct;

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
    temp = load(hypfile);
    clear temp;
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

    % true data
    posterior = mean(ratesPosterior.(name),1);%temporal mean
    pomean = mean(posterior);
    % urdme data
    posterior_urdme = mean(ratesPosterior_urdme.(name),1);
    pomean_urdme = mean(posterior_urdme);
    switch name
      case {'F0' 'HOSP' 'IC_HOSP' 'IFR' 'E2I'}
        prior = prior.*100;
        prmean = prmean*100;
        posterior = posterior.*100;
        pomean = pomean*100;
        posterior_urdme = posterior_urdme.*100;
        pomean_urdme = pomean_urdme*100;
      case {'sigma' 'gammaI' 'gammaH' 'gammaW'}
        prior = 1./prior;
        prmean = 1/prmean;
        posterior = 1./posterior;
        pomean = mean(posterior); % (need not commute)
        posterior_urdme = 1./posterior_urdme;
        pomean_urdme = mean(posterior_urdme); % (need not commute)
      case 'half_life'
        prior = prior*24;
        prmean = prmean*24;
        posterior = posterior*24;
        pomean = pomean*24;
        posterior_urdme = posterior_urdme*24;
        pomean_urdme = pomean_urdme*24;
    end

    % posterior | real data
    [fp,xip] = ksdensity(posterior,'Function','pdf');
    xx = [xip,flip(xip,1)];
    ff = [fp,flip(fp,1)];
    h_post = fill(xx,ff,[0 0 1]);
    hold on
    set(h_post,'facealpha',.1)

    % posterior | URDME data
    [fp,xip] = ksdensity(posterior_urdme,'Function','pdf');
    xx = [xip,flip(xip,1)];
    ff = [fp,flip(fp,1)];
    h_boot = fill(xx,ff,[0 1 0]);
    hold on
    set(h_boot,'facealpha',.1)

    % prior
    [f,xi] = ksdensity(prior,'Function','pdf');
    xx = [xi flip(xi,1)];
    ff = [f flip(f,1)];
    h_prior = fill(xx,ff,[1 0 0]);
    set(h_prior,'facealpha',.1)

    xlabel(rateNamesSubset_tex{i},'interpreter','latex');

    % mean vertical line, with annotation.
    hline_prior = xline(prmean,'--','LineWidth',1,'Color', [1 0.2 0.2]);%[0.9844 0.5820 0.5703]);
    hline_post = xline(pomean,'-','LineWidth',2,'Color', [0 0.4470 0.7410]);
    str = ['$=$ ' num2str(round(pomean,2,'significant'))];
    b(i) = annotation('textbox','String',str,'Position',ha(i).Position,...
                      'Vert','top','HorizontalAlignment','right',...
                      'FitBoxToText','on','interpreter','latex',...
                      'BackgroundColor',[1 1 1], ...
                      'FaceAlpha',0.01,...
                      'EdgeColor','None','Margin',0,'FontSize',8,...
                      'Color',[0 0.4470 0.7410]);

    hline_boot = xline(pomean_urdme,':b','LineWidth',2,'Color', [0.4102 0.6992 0.4844]);
    str_urdme = ['$=$ ' num2str(round(pomean_urdme,2,'significant'))];
    b(i) = annotation('textbox','String',str_urdme,'Position',ha(i).Position,...
                      'Vert','mid','HorizontalAlignment','right',...
                      'FitBoxToText','on','interpreter','latex',...
                      'BackgroundColor',[1 1 1], ...
                      'FaceAlpha',0.01,...
                      'EdgeColor','None','Margin',0,'FontSize',8,...
                      'Color',[0.4102 0.6992 0.4844]);



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
      case 'HOSP', xlim([0 10]);
      case 'IC_HOSP', xlim([0 40]);
    end
    set(gca,'ytick',[]);
    set(get(gca,'Xruler'),'Exponent',0)
    hold off


end

legend([h_post(1), h_prior(1), h_boot(1)], ...
       'Posterior','Prior','Bootstrapped Posterior', ...
       'Orientation','horizontal','interpreter','latex', ...
       'Position',[0.47 0.95 0.1 0.03],'FontSize',14);


% legend('Posterior','Bootstrapped Posterior', 'Prior', ...
%   'Orientation','horizontal','interpreter','latex', ...
%   'Position',[0.47 0.95 0.1 0.03],'FontSize',14);
legend boxoff




%% publication-friendly formats
abspath = mfilename('fullpath');
figname = [abspath(1:end-29) '/posterior_urdmesweden.pdf'];
if savetofile
    % figure(1),
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'Position',[100 100 700 400]);
    print('-dpdf', figname)
    if verb
        disp(['saved figure: ' figname])
    end
else
    if verb
        disp(['didn''t save figure: ' figname])
    end
end