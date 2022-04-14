% DEADSPLIT generates CrI and Figures concerning the CFR and D source
% compartment split.
%
% WARNING in order to generate the data you must first checkout the H2
% branch and generate the data file. Then preferably pull it into the
% current branch for furhter processing by this script. The reason for
% this procedure is that H2 introduces breaking changes to the
% 'legacy' solutions by adding new compartments. A full merge can be
% done, but is not a priority as of writing.
%
% Robin Marin 22-03-22

% Load simulation
ending = '1';
try
  load(['weekly/save/runs/slam210531_Sweden_monthly_' ending])
catch
    error('you tried to run the script without the necessary simulation file');
end

spectral; % fetch ODE CFR:s and D_SS_IHW

% load data
Data = loadData('C19');
Data = polishData(Data,'D','Dinc',1);
Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});

% compare to what date?
stopdate=210323;

% Get the correct fields
% 1:H 2:W 3:H2 4:D_I 5:D_H 6:D_W 7:PHI 8:I 9:R 10:Ic 11:Hc 12:Wc
base = size(Z,1) - 12;
sel_H = base + [1 3];
sel_D = base + [4 5 6];
sel_Ic = base + [10];
sel_Hc = base + [11];
sel_Wc = base + [12];

%% Does H+H2 match H_data?
%
figure(1), clf

sel_H = base + [1 3];
h1=plot(squeeze(sum(Z(sel_H,:,:),1)),'k');
hold on
h2=plot(squeeze(sum(Z(sel_H(1),:,:),1)),'b');
h3=plot(squeeze(sum(Z(sel_H(2),:,:),1)),'g');
h4=plot(sum(Data.H(tspan_data,:),2),'--r');
legend([h1(1) h2(1) h3(1) h4(1)], ...
  {'H+H2' 'H' 'H2' 'data'},...
'Location','northwest')
title('H_{data} vs H + H_2')

%% Does D_I+D_H+D_W match D_data?
figure(2), clf;
h1 = plot(squeeze(sum(Z(sel_D,:,:),1)),'k');
hold on
h2 = plot(squeeze(sum(Z(sel_D(1),:,:),1)),'b');
h3 = plot(squeeze(sum(Z(sel_D(2),:,:),1)),'y');
h4 = plot(squeeze(sum(Z(sel_D(3),:,:),1)),'g');
h5 = plot(sum(Data.D(tspan_data,:),2),'--r');
xh = xline(find(Data.date(tspan_data) == stopdate),':');
legend([h1(1) h2(1) h3(1) h4(1) h5(1) xh(1)], ...
  {'sum(D_*)' 'D_I' 'D_H' 'D_W' 'data' '210323'},...
'Location','northwest')
title('D_{data} vs D_I + D_H + D_W')




%% Get the distribution of {D_I ...} at date of comparison
D_post_IHW=Z(sel_D,Data.date(tspan_data) == stopdate,:);
D_post_IHW_tot = sum(D_post_IHW);
D_post_IHW_norm = 100*(D_post_IHW./D_post_IHW_tot);

%%
D_post_IHW_sig = squeeze(covZ.stdZ(sel_D,Data.date(tspan_data) == stopdate,:));
D_post_IHW_mu = squeeze(D_post_IHW);

xx = linspace(0,1e4,1e3);
[mu_IHW, sig_IHW] = l_MoG(D_post_IHW,D_post_IHW_sig)

pd = normpdf(xx,mu_IHW,sig_IHW);


% quantiles on data
D_post_IHW_q = quantile(squeeze(D_post_IHW),[0.05/2 0.32/2 0.5 1-0.32/2 1-0.05/2],2)
% and normalized
D_post_IHW_norm_q = quantile(squeeze(D_post_IHW_norm),[0.05/2 0.32/2 0.5 1-0.32/2 1-0.05/2],2)


colors = {'#0072BD' '#D95319' '#EDB120'};
figure(3), clf, hold on
for k = 1:3
  histogram(D_post_IHW(k,:),10,'Normalization','pdf', 'Facecolor',colors{k},'Facealpha',0.2)
  plot(xx,pd(k,:),'Color',colors{k},'LineWidth',3,'HandleVisibility','off')
  xline(D_SS_IHW(k),'--','Color',colors{k},'HandleVisibility','off','LineWidth',3)
  xline(mu_IHW(k)-2*sig_IHW(k),'-','Color',colors{k},'HandleVisibility','off','LineWidth',2)
  xline(mu_IHW(k)+2*sig_IHW(k),'-','Color',colors{k},'HandleVisibility','off','LineWidth',2)

end

hold off
legend({'D_I' 'D_H' 'D_W'})
xlim([1e1 1e4])
xlabel('Total # of dead per compartment source')
title('Density at date: D_I + D_H + D_W')

%% X --> D (eventually) = X/D
figure(4), clf
tiledlayout(1,3)


% I -> D
nexttile;
I_post_sig = covZ.stdZ(sel_Ic,Data.date(tspan_data) == stopdate,:);
D_post_sig = covZ.stdZ(sel_D,Data.date(tspan_data) == stopdate,:);
I_post = Z(sel_Ic,Data.date(tspan_data) == stopdate,:);

histogram(100*squeeze(sum(D_post_IHW([1 2 3],:,:),1)./I_post),'normalization','pdf'), hold on

% Z = X/Y, where X~Normal and Y~normal
% is a weird distribution, see https://rstudio-pubs-static.s3.amazonaws.com/287838_7c982110ffe44d1eb5184739c5724926.html
% but is ~Normal if delta_Y := CV(Y) < 0.1
[mu_Id, sig_Id] = l_MoG(sum(D_post_IHW([1 2 3],:,:),1),sqrt(sum(D_post_sig([1 2 3],:,:).^2,1)), I_post, I_post_sig);

mu_Id = 100*mu_Id; sig_Id = 100*sig_Id; % scale for [%] in plot


xx=100*linspace(0,0.03,100);
plot(xx,normpdf(xx,mu_Id,sig_Id))

% 95% CrI for MoG
xline(1*mu_Id,'-','Linewidth',2,'Color',colors{2})
xline(1*(mu_Id-2*sig_Id),'--r','Linewidth',2,'Color',colors{2})
xline(1*(mu_Id+2*sig_Id),'--r','Linewidth',2,'Color',colors{2})

% 95% CrI for ODE
xline(ID_q(3),'-','Linewidth',2)
xline(ID_q(1),'--','Linewidth',2)
xline(ID_q(5),'--','Linewidth',2)
title('I \rightarrow D (eventually)')
xlabel('%')



% H -> D
nexttile;
H_post = Z(sel_Hc,Data.date(tspan_data) == stopdate,:);
H_post_sig = covZ.stdZ(sel_Hc,Data.date(tspan_data) == stopdate,:);
[mu_Hd, sig_Hd] = l_MoG(sum(D_post_IHW([2 3],:,:),1),sqrt(sum(D_post_sig([2 3],:,:).^2,1)),...
                        H_post, H_post_sig);
mu_Hd = 100*mu_Hd; sig_Hd = 100*sig_Hd;
histogram(100*squeeze(sum(D_post_IHW([2 3],:,:),1)./H_post),'normalization','pdf'), hold on

xx=100*linspace(0,0.4,100);
plot(1*xx,normpdf(xx,mu_Hd,sig_Hd))

% 95% CrI for MoG
xline(1*mu_Hd,'-','Linewidth',2,'Color',colors{2})
xline(1*(mu_Hd-2*sig_Hd),'--r','Linewidth',2,'Color',colors{2})
xline(1*(mu_Hd+2*sig_Hd),'--r','Linewidth',2,'Color',colors{2})

% 95% CrI for ODE
xline(HD_q(3),'-','Linewidth',2)
xline(HD_q(1),'--','Linewidth',2)
xline(HD_q(5),'--','Linewidth',2)
title('H \rightarrow D (eventually)')
xlabel('%')


% W -> D
nexttile;
W_post = Z(sel_Wc,Data.date(tspan_data) == stopdate,:);
W_post_sig = covZ.stdZ(sel_Wc,Data.date(tspan_data) == stopdate,:);
[mu_Wd, sig_Wd] = l_MoG(sum(D_post_IHW([3],:,:),1),sqrt(sum(D_post_sig([3],:,:).^2,1)),...
                        W_post, W_post_sig);
mu_Wd = 100*mu_Wd; sig_Wd = 100*sig_Wd;


histogram(100*squeeze(sum(D_post_IHW(3,:,:),1)./W_post),'normalization','pdf'), hold on


xx=100*linspace(0.1,0.9,100);
plot(xx,normpdf(xx,mu_Wd,sig_Wd))

% 95% CrI for MoG
xline(1*mu_Wd,'-','Linewidth',2,'Color',colors{2})
xline(1*(mu_Wd-2*sig_Wd),'--r','Linewidth',2,'Color',colors{2})
xline(1*(mu_Wd+2*sig_Wd),'--r','Linewidth',2,'Color',colors{2})

% 95% CrI for ODE
xline(WD_q(3),'-','Linewidth',2)
xline(WD_q(1),'--','Linewidth',2)
xline(WD_q(5),'--','Linewidth',2)
title('W \rightarrow D (eventually)')
%legend({'P. pred' 'T. mean' 'T. 95% CrI'})
xlabel('%')

% %%
% figure(5), clf,
%
% histogram(100*squeeze(sum(D_post_IHW(3,:,:),1)./W_post),'normalization','pdf'), hold on
%
%
% xx=100*linspace(0.1,0.9,100);
% plot(xx,normpdf(xx,mu_Wd,sig_Wd))
%
% % 95% CrI for MoG
% xline(1*mu_Wd,'-','Linewidth',2,'Color',colors{2})
% xline(1*(mu_Wd-2*sig_Wd),'--r','Linewidth',2,'Color',colors{2})
% xline(1*(mu_Wd+2*sig_Wd),'--r','Linewidth',2,'Color',colors{2})
%
% % 95% CrI for ODE
% xline(WD_q(3),'-','Linewidth',2)
% xline(WD_q(1),'--','Linewidth',2)
% xline(WD_q(5),'--','Linewidth',2)
% title('W \rightarrow D (eventually) [ZOOM]')
% %legend({'P. pred' 'T. mean' 'T. 95% CrI'})
% xlabel('%')
% xlim([30,60])
% hold off

%% Format intervals for LaTeX use
disp('D split distribution')
D_Kalman_IHW = l_formatCI(mu_IHW + [-2*sig_IHW, - sig_IHW, zeros(3,1), sig_IHW, 2*sig_IHW], 4)
D_SS_IHW

disp('CFR calculation')
CFR_Kalman_IHW = l_formatCI([mu_Id;mu_Hd;mu_Wd] + ...
                          [-2*sig_Id, - sig_Id, 0, sig_Id, 2*sig_Id;...
                           -2*sig_Hd, - sig_Hd, 0, sig_Hd, 2*sig_Hd;...
                           -2*sig_Wd, - sig_Wd, 0, sig_Wd, 2*sig_Wd],2)
ID_ODE = l_formatCI([ID_q;HD_q;WD_q],2)



return;

%%
function [mu,sig] = l_MoG(X_mu,X_sig,Y_mu,Y_sig)
%L_MOG computes the Mixture of Gaussian approximation. For
% case 1 X~Norm
% or
% case 2: Z=X/Y, where X~Norm, Y~Norm, and approximate Z~Norm

if nargin == 4 % MoG approximation of a ~ Cauchi
  mu_Z = X_mu ./ Y_mu;
  delta_Y = Y_sig./Y_mu; % CoV
  rho_Z = Y_sig ./ X_sig;

  sig_Z = sqrt(delta_Y.^2 .*( rho_Z.^(-2) + mu_Z));

  % format to proper shape
  mu_Z = squeeze(mu_Z)';
  sig_Z = squeeze(sig_Z)';
else
  mu_Z = squeeze(X_mu);
  sig_Z = squeeze(X_sig);
end

mu = mean(mu_Z,2);
getsig = @(mu,mu_,sig_) sqrt(mean(mu_.^2 + sig_.^2,2) - mu.^2);
sig = getsig(mu,mu_Z,sig_Z);

end

function [CI] = l_formatCI(X,Nsig)
%L_FORMATCI converts numeric array X into a string \CI{...} format ready
%for LaTeX use.

CI = cell(size(X,1),1);
X = round(X,Nsig,'significant');
for k = 1:size(X,1)
 CI{k} = ['\text{\CIedge' ...
   '{' num2str(X(k,1)) '}'...
   '{' num2str(X(k,2)) '}'...
   '{' num2str(X(k,3)) '}'...
   '{' num2str(X(k,4)) '}'...
   '{' num2str(X(k,5)) '}}'];

end
end
