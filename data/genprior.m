%genPrior Generate prior.
%   This script generates the constants required for the EngEr
%   Covid-19 model prior.
%
% References:
%
% [1] https://www.socialstyrelsen.se/globalassets/1-globalt/
%     covid-19-statistik/vard-och-covid-19/statistik-covid19-vardforlopp.xlsx,
%     retrieved 21-05-07
% [2] W. Dhouib et al.: "The incubation period during the pandemic of
%     COVID-19: a systematic review and meta-analysis", Syst. Rev.10(101), 2021
%     https://systematicreviewsjournal.biomedcentral.com/articles/10.1186/s13643-021-01648-y
% [3] M. Alene et al.: "Magnitude of asymptomatic COVID-19 cases
%     throughout the course of infection: A systematic review and
%     meta-analysis", PloS ONE, 2021
%     https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0249090
% [4] https://portal.icuregswe.org/siri/en/report/corona.mort?filter=d5d1d908-cd6a-8db0-e5ec-8c76b834d200
%     retrieved 21-05-03

% S. Engblom 2021-05-07 (inspired by generateFcases.m)

% illustrate or don't
illustrate = true;

% (A) rates of the model

% rates sigma, gammaI, gammaH, gammaW
[sigmaPQ,sigmaHyp,gammaIPQ,gammaIHyp] = l_sigma_gammaI(illustrate);
[gammaHPQ,gammaHHyp,gammaWPQ,gammaWHyp] = l_gammaHW(illustrate);

% means are useful too:
sigmamean = l_scaledbetastats(sigmaHyp,sigmaPQ);
gammaImean = l_scaledbetastats(gammaIHyp,gammaIPQ);
gammaHmean = l_scaledbetastats(gammaHHyp,gammaHPQ);
gammaWmean = l_scaledbetastats(gammaWHyp,gammaWPQ);

% (B) fractions of the model

% fraction E --> I here taken from [3]: 95% CI: 0.75 [0.62 0.84],
% determine the Beta distribution that fits this:
fun = @(abc)abc(3)+(1-abc(3))* ...
      [betastat(abc(1),abc(2)) betacdf([0.62 0.84],abc(1),abc(2))]- ...
      [75 2.5 97.5]*1e-2;
if illustrate
  ABC = fsolve(fun,[3 2 0.01]);
  figure(8), clf,
  xx = linspace(ABC(3),1);
  plot(xx,l_scaledbetapdf(xx,ABC(1:2),[ABC(3) 1]));
  title('E2I: fraction E $\to$ I','interpreter','latex');
else
  ABC = fsolve(fun,[3 2 0.01],optimoptions(@fsolve,'Display','off'));
end
E2IPQ = [ABC(3) 1]; % support on (c,1)
E2IHyp(1) = ABC(1);
E2IHyp(2) = ABC(2);
E2Imean = l_scaledbetastats(E2IHyp,E2IPQ);
% (F0 = E2I)

% fraction E --> D (the Infection Fatality Rate, IFR)
IFRPQ = [0 0.02]; % supp = [0,2]%
IFRHyp = [2 4]; % mean = 2/3%
IFRmean = l_scaledbetastats(IFRHyp,IFRPQ);
% (IFR is not used directly but later determines the rate I_MORT: I
% --> D after taking {H,W} --> D into account)

% k (scale factor between phi_sew and phi s.t. phi_sew = k_sew*phi

%load k_sew into workspace, see sew_plots
abspath = mfilename('fullpath');
load([abspath(1:end-13) 'sewage/k_sew_timeseries'])

% Set the desired mean and variance for beta prior
mu_k_sew = 2.5*MK_mean;
sigma2_k_sew = (4*MK_std)^2;

%% This would take the mean and 2*variance from a uniform
% distribution between the min and max of k_sew. However, does not
% result in a "bubbly" beta distribution
%a = 0;
%b = 0.015;
%mu_k_sew = (a+b)/2;
%sigma2_k_sew = 2*(a-b)^2/12;

% Set limits of beta prior
k_sew_PQ = [0 0.015];

% Calculate hyperparameters for beta prior that give the desired
% mean and variance specified above:
lambda_k_sew = (mu_k_sew - k_sew_PQ(1))*(k_sew_PQ(2)-mu_k_sew)/sigma2_k_sew-1;

alpha_k_sew = lambda_k_sew*(mu_k_sew - k_sew_PQ(1))/(k_sew_PQ(2) - k_sew_PQ(1));
beta_k_sew = lambda_k_sew*(k_sew_PQ(2)-mu_k_sew)/(k_sew_PQ(2)-k_sew_PQ(1)); 

k_sew_Hyp = [alpha_k_sew beta_k_sew];

if illustrate
  figure(12), clf
  XXX = 0:0.00001:0.015;
  YYY = l_scaledbetapdf(XXX, k_sew_Hyp, k_sew_PQ);
  plot(XXX, YYY, 'Color','red')
  hold on
  histogram(k_sew,20, 'Normalization', 'pdf', 'FaceColor','b')
  legend('proposed beta-fit', '"data"')
  title('$k_{sewage}$','interpreter','latex');
end
% (Should be same as mu_k_sew above)
k_sew_mean = l_scaledbetastats(k_sew_Hyp, k_sew_PQ);

% SIR_MORT: mortality under IC (W --> D)
SIR = loadData('SIR_swe_mort'); % (source: [4])
SIR_MORT = mean(SIR.SIR_MORT);
% (F4 = SIR_MORT)

% HOSP_MORT: mortality under hospitalization (H --> D)
hosp_mort = 5729/46236; % (source: [1])
% (note: relative to the IC mortality)
HOSP_MORT = mean(hosp_mort./SIR.SIR_MORT);
% (F3d = SIR_MORT*HOSP_MORT)

% IC_HOSP: fraction H --> W; this fraction is very elusive so here we
% make a coarse fit based on aggregated data
Data = loadData('C19');
Data = polishData(Data,{},{},1);
Data = smoothData(Data,{'H' 'W'},{[] []});

% per day and region, but remove small values of W
HWfac = reshape(Data.H./Data.W,[],1);
HWfac(Data.W < 0.5*nanmean(Data.W(:))) = NaN;
HWfac(isnan(HWfac)) = [];
IC_HOSP_ = gammaWmean./(gammaHmean*HWfac);
ABC = quantile(IC_HOSP_,[0.005 0.5 0.995]);

% Note, there exists an alternative solution for IC_HOSP which is
% achived by the following steps: per day, but all of Sweden
% HWfac = sum(Data.H,2)./sum(Data.W,2);
% IC_HOSP_ = gammaWmean./(gammaHmean*HWfac);
% ABC = quantile(IC_HOSP_,[0.05 0.5 0.95]);
% Beta with support at the 90% CI and with mean at the median
% IC_HOSPPQ = ABC([1 3]);
% IC_HOSPHyp(1) = 2;
% IC_HOSPHyp(2) = IC_HOSPHyp(1)/((ABC(2)-IC_HOSPPQ(1))/diff(IC_HOSPPQ))- ...
%    IC_HOSPHyp(1);
% IC_HOSPmean = l_scaledbetastats(IC_HOSPHyp,IC_HOSPPQ);





% Beta-fit with support at the 99% CI and with mean at the median
IC_HOSPPQ = ABC([1 3]);
IC_HOSPHyp(1) = 2;
IC_HOSPHyp(2) = IC_HOSPHyp(1)/((ABC(2)-IC_HOSPPQ(1))/diff(IC_HOSPPQ))- ...
    IC_HOSPHyp(1);
IC_HOSPmean = l_scaledbetastats(IC_HOSPHyp,IC_HOSPPQ);

% The above expression comes from a quasi steady-state assumption:
% suppose that H ~ HWfac*W when W is stationary. This means that
%   W_in = H*IC_HOSP*gammaH ~ W*gammaW = W_out
% such that
%   IC_HOSP ~ W*gammaW/(gammaH*H) ~ gammaW/(gammaH*HWfac).

% (F3 = IC_HOSP)

if illustrate
  figure(10), clf,
  xx = linspace(IC_HOSPPQ(1),IC_HOSPPQ(2));
  histogram(IC_HOSP_,100,'Normalization','pdf'); hold on
  plot(xx,l_scaledbetapdf(xx,IC_HOSPHyp,IC_HOSPPQ));
  legend('Sample','Proposed beta-fit');
  title('IC\_HOSP: H $\to$ W','interpreter','latex');
end

% At this point we have the relation
%   IFR = E2I*I_MORT+[H]*HOSP_MORT*SIR_MORT+[W]*SIR_MORT.
% We find the effective values of [H] and[W] by recognizing a
% geometric series,
%   [H] = E2I*HOSP*(1+x+x^2+...) with x = IC_HOSP*(1-SIR_MORT),
%   [H] = E2I*HOSP/(1-x),
% and
%   [W] = E2I*HOSP*IC_HOSP*(1+x+x^2+...),
%   [W] = E2I*HOSP*IC_HOSP/(1-x).
% To sum up,
%   IFR = E2I*[I_MORT+(HOSP_MORT+IC_HOSP)*SIR_MORT*HOSP/(1-x)].

% HOSP: fraction I --> H
IFR_ = l_scaledbetarnd(IFRHyp,IFRPQ,1,1e5);
E2I_ = l_scaledbetarnd(E2IHyp,E2IPQ,1,1e5);
IC_HOSP_ = l_scaledbetarnd(IC_HOSPHyp,IC_HOSPPQ,1,1e5);
x = IC_HOSP_*(1-SIR_MORT);
% Connecting model: total mortality in I can be written as
% I_MORT+H_MORT+W_MORT with H_MORT = HOSP_MORT*SIR_MORT*HOSP/(1-x) and
% W_MORT = IC_HOSP*SIR_MORT*HOSP/(1-x). We assume there is a relation
% of the form I_MORT = I2HWfac*(H_MORT+W_MORT) such that the total
% mortality of I becomes (1+I2HWfac)*(H_MORT+W_MORT).
%
% Now, I2HWfac can be estimated:
% Hmort_tot = 5729; % (source: [1], date 210504)
% Wmort_tot = 1457;
% Data = loadData('FHM_swe');
% ij = find(Data.date == 210504);
% Imort_tot = nansum(Data.D(1:ij))-Hmort_tot-Wmort_tot;
% I2HWfac = Imort_tot/(Hmort_tot+Wmort_tot); % ~0.9763

% simulate the resulting model and fit a beta:
I2HWfac_ = 2*betarnd(2,2,1,1e5); % supp(I2HWfac) = (0,2)
HOSP_ = (IFR_./(E2I_.*(1+I2HWfac_)))./ ...
        (((HOSP_MORT+IC_HOSP_)*SIR_MORT).*(1-x));
HOSPPQ = [0 max(HOSP_)*(1+1e-6)];
HOSPHyp = betafit(HOSP_/HOSPPQ(2));
HOSPmean = l_scaledbetastats(HOSPHyp,HOSPPQ);

if illustrate
  figure(9), clf,
  HOSPPQold = [0 (IFRmean/E2Imean)/ ...
               ((HOSP_MORT+IC_HOSPmean)*SIR_MORT)* ...
               (1-IC_HOSPmean*(1-SIR_MORT))];
  HOSPHypold = [2 2];
  HOSPmeanold = l_scaledbetastats(HOSPHypold,HOSPPQold);
  xx = linspace(0,max(HOSPPQ(2),HOSPPQold(2)));
  histogram(HOSP_,40,'Normalization','pdf'); hold on
  plot(xx,l_scaledbetapdf(xx,HOSPHyp,HOSPPQ));
  plot(xx,l_scaledbetapdf(xx,HOSPHypold,HOSPPQold));
  legend('Monte Carlo','Proposed beta-fit','Current prior');
  title('HOSP: fraction I $\to$ H','interpreter','latex');
end
% (F2 = HOSP)

% Through actual sampling of the prior, we may subsequently deduce a
% suitable value of I_MORT which agrees with the sampled value for the
% IFR. That is,
%   I_MORT = IFR/E2I-(HOSP_MORT+IC_HOSP)*SIR_MORT*HOSP/(1-x),
% with x = IC_HOSP*(1-SIR_MORT).
% (F2d = I_MORT)

% for convenience:
HOSP_ = l_scaledbetarnd(HOSPHyp,HOSPPQ,1,1e5);
I_MORT_ = max(IFR_./E2I_-(HOSP_MORT+IC_HOSP_)*SIR_MORT.*HOSP_./(1-x),0);
I_MORTmean = mean(I_MORT_);
if illustrate
  figure(11), clf,
  histogram(I_MORT_,80,'Normalization','pdf');
  ylim([0 200]);
  legend('Monte Carlo');
  title('I\_MORT: I $\to$ D','interpreter','latex');
end

% pack up as an output struct
c19prior.sigmaPQ = sigmaPQ;
c19prior.sigmaHyp = sigmaHyp;
c19prior.sigmamean = sigmamean;
c19prior.gammaIPQ = gammaIPQ;
c19prior.gammaIHyp = gammaIHyp;
c19prior.gammaImean = gammaImean;

c19prior.gammaHPQ = gammaHPQ;
c19prior.gammaHHyp = gammaHHyp;
c19prior.gammaHmean = gammaHmean;
c19prior.gammaWPQ = gammaWPQ;
c19prior.gammaWHyp = gammaWHyp;
c19prior.gammaWmean = gammaWmean;

c19prior.E2IPQ = E2IPQ;
c19prior.E2IHyp = E2IHyp;
c19prior.E2Imean = E2Imean;
c19prior.IFRPQ = IFRPQ;
c19prior.IFRHyp = IFRHyp;
c19prior.IFRmean = IFRmean;

c19prior.k_sew_PQ = k_sew_PQ;
c19prior.k_sew_Hyp = k_sew_Hyp;
c19prior.k_sew_mean = k_sew_mean;

% sole two point estimates:
c19prior.SIR_MORT = SIR_MORT;
c19prior.HOSP_MORT = HOSP_MORT;

c19prior.IC_HOSPPQ = IC_HOSPPQ;
c19prior.IC_HOSPHyp = IC_HOSPHyp;
c19prior.IC_HOSPmean = IC_HOSPmean;
c19prior.HOSPPQ = HOSPPQ;
c19prior.HOSPHyp = HOSPHyp;
c19prior.HOSPmean = HOSPmean;

% sampled later on the fly, but for convenience:
c19prior.I_MORTmean = I_MORTmean;

% add date of calculation
c19prior.rev = date;
% END OF SCRIPT

% ----------------------------------------------------------------------
function p = l_scaledbetapdf(x,AB,LU)
%L_SCALEDBETAPDF Density of scaled & shifted beta density.
%   P = L_SCALEDBETAPDF(X,AB,LU) computes the density P for points X
%   of the beta density of parameters AB(1) and AB(2). The density has
%   support in (LU(1),LU(2)).
p = betapdf((x-LU(1))/diff(LU),AB(1),AB(2))/diff(LU);
end
% ----------------------------------------------------------------------
function r = l_scaledbetarnd(AB,LU,M,N)
%L_SCALEDBETARND Sample from scaled & shifted beta density.
r = LU(1)+diff(LU)*betarnd(AB(1),AB(2),M,N);
end
% ----------------------------------------------------------------------
function [m,v] = l_scaledbetastats(AB,LU)
%L_SCALEDBETASTATS Mean and variance of scaled & shifted beta density.
[m,v] = betastat(AB(1),AB(2));
m = LU(1)+diff(LU)*m;
v = diff(LU)^2*v;
end
% ----------------------------------------------------------------------
function [sigmaPQ,sigmaHyp,gammaIPQ,gammaIHyp] = l_sigma_gammaI(ill)
%L_SIGMA_GAMMAI Suitable priors for sigma and gammaI/gammaA.

% mean incubation time is 6.2 with [5.4,7.0] a 95% CI (source: [2])
sigmaPQ = 1./[7.0 5.4]; % prior support
sigmaMU = 1/6.2; % prior mean

sigmaHyp(1) = 2;
sigmaHyp(2) = sigmaHyp(1)/((sigmaMU-sigmaPQ(1))/diff(sigmaPQ))-sigmaHyp(1);

if ill
  figure(6), clf, hold on
  sample = l_scaledbetarnd(sigmaHyp,sigmaPQ,1,1e5);
  sample = exprnd(1./sample);
  days = 0:25;
  sampledens = histc(sample,[days inf]);
  sampledens_sigma = sampledens(1:end-1)/sum(sampledens);

  xline(quantile(sample,0.95),'b--','Upper 95\% CI [$\sigma$]', ...
        'HandleVisibility','off','interpreter','latex');
end

% less information is available on gammaI
gammaIPQ = 1./[10 4]; % prior support
gammaIMU = 1/7; % prior mean

gammaIHyp(1) = 2;
gammaIHyp(2) = gammaIHyp(1)/((gammaIMU-gammaIPQ(1))/diff(gammaIPQ))-gammaIHyp(1);

if ill
  sample = l_scaledbetarnd(gammaIHyp,gammaIPQ,1,1e5);
  sample = exprnd(1./sample);
  sampledens = histc(sample,[days inf]);
  sampledens_gammaI = sampledens(1:end-1)/sum(sampledens);

  bar(days,[sampledens_sigma; sampledens_gammaI]');
  legend('$\sigma$','$\gamma_I$','Location','N','interpreter','latex');
  xlabel('days');

  xline(quantile(sample,0.95),'b--','Upper 95\% CI [$\gamma_I$]', ...
        'HandleVisibility','off','interpreter','latex');
  title('Exit time: incubation time and symptomatic time');
  hold off

  figure(7), clf, hold on
  xx = linspace(1/10,1/4);
  plot(1./xx,1./xx.^2.*l_scaledbetapdf(xx,sigmaHyp,sigmaPQ),'b--');
  plot(1./xx,1./xx.^2.*l_scaledbetapdf(xx,gammaIHyp,gammaIPQ),'r--');
  xlabel('days');
  title('Hyperparameter density: incubation and symptomatic');
  hold off
end

end
% ----------------------------------------------------------------------
function [gammaHPQ,gammaHHyp,gammaWPQ,gammaWHyp] = l_gammaHW(ill)
%L_GAMMAHW Suitable priors for gammaH and gammaW from data.

% (1) gammaH: empirical fit (source: [1])
NobsH = 40507;
daysHfrac = [2.6 10.0 11.3 10.1 8.8 7.7 6.6 6.0 4.9 4.2 3.5 ...
             3.1 2.6 2.3 2.1 1.7 1.3 1.2 1.0 0.8 0.9 0.8 0.6 ...
             0.5 0.5 0.5 0.4 0.3 0.3 0.3 3.0];
daysHfrac = daysHfrac/sum(daysHfrac);
daysH = 1:31; % (note: final fraction is > 30 days)

% lognormal fit (only for visualization)
[phat,phatci] = lognfit(daysH,0.05,[zeros(1,numel(daysH)-1) 1],daysHfrac);
[gammaH,vargammaH] = lognstat(phat(1),phat(2));
gammaH = 1/gammaH; % = 1/8.99

% Bayesian inference via gamma conjugate prior
prior = [1 1/7]; % [alpha beta] mean = 1/7 days, std = 1/7 days
% (here one should normally use the mean, but this does not allow for
% right-censoring, hence we use the mean from an ML-fit instead)
lam = expfit(daysH,0.05,[zeros(1,numel(daysH)-1) 1],daysHfrac);
posteriorH = [prior(1)+NobsH 1/(prior(2)+NobsH*lam)];
% (formula from https://en.wikipedia.org/wiki/Exponential_distribution#Bayesian_inference)

% Bayesian prior (incl simple beta-fit) vs current prior
if ill, figure(1), clf, hold on, end
xx = gaminv([5e-2/2 1-5e-2/2],posteriorH(1),posteriorH(2));
xx = linspace(xx(1),xx(2));
if ill
  plot(1./xx,1./xx.^2.*gampdf(xx,posteriorH(1),posteriorH(2)),'b');
end
% here is the proposed new fit:
m = mean(xx([1 end]));
s = 0.25*diff(xx([1 end]));
xx = linspace(m-4*s,m+4*s);
gammaHPQ = xx([1 end]);
gammaHHyp = [3 3];
if ill
  plot(1./xx,1./xx.^2.*l_scaledbetapdf(xx,gammaHHyp,gammaHPQ),'k--');
end

% current prior
gammaHPQold = [1/9 1/6];
gammaHMU = 1/8.26;
gammaHHypold(1) = 2;
gammaHHypold(2) = gammaHHypold(1)/((gammaHMU-gammaHPQold(1))/ ...
                             diff(gammaHPQold))-gammaHHypold(1);
xx = linspace(gammaHPQold(1),gammaHPQold(2));
if ill
  plot(1./xx,1./xx.^2.*l_scaledbetapdf(xx,gammaHHypold,gammaHPQold),'r')
  title('Hyperparameter density: hospital');
  legend('Bayesian $\Gamma$-fit 95\% CI','Proposed beta-fit', ...
         'Current prior','Location','NW','interpreter','latex');
  xlabel('days');
end

if ill
  % illustration of fit
  figure(2), clf
  % original gamma-fit
% $$$ sample = exprnd(1./gamrnd(posteriorH(1),posteriorH(2),1,1e5));
% $$$ sampledens = histc(sample,[0 daysH(1:end-1) inf]);
% $$$ sampledens = sampledens(1:end-1)/sum(sampledens);

  % relaxed beta-fit
  sample = l_scaledbetarnd(gammaHHyp,gammaHPQ,1,1e5);
  sample = exprnd(1./sample);
  sampledens = histc(sample,[0 daysH(1:end-1) inf]);
  sampledens = sampledens(1:end-1)/sum(sampledens);

  % current prior
  sample_ = gammaHPQold(1)+diff(gammaHPQold)* ...
            betarnd(gammaHHypold(1),gammaHHypold(2),1,1e5);
  sample_ = exprnd(1./sample_);
  sampledens_ = histc(sample_,[0 daysH(1:end-1) inf]);
  sampledens_ = sampledens_(1:end-1)/sum(sampledens_);

  bar(daysH,[daysHfrac; sampledens; sampledens_]'), hold on
  plot(daysH,lognpdf(daysH,phat(1),phat(2)),'ro-');
  title('Exit time hospital');
  legend('Data','Proposed fit','Current prior','lognormal');
  xlabel('days');
end

% (2) gammaW: as above...

% (2.1) ...but firstly, data is "gammaH+gammaW+gammaH" so we need to
% subtract gammaH+gammaH (understood as independent waiting times)

% impute the empirical distribution to the right
daysH_ = 31:60;
daysHfrac_ = [daysHfrac(1:end-1) lognpdf(daysH_,phat(1),phat(2))];
daysHfrac_ = [1-sum(daysHfrac_) daysHfrac_];
daysH_ = [0 daysH(1:end-1) daysH_]; % note: add support at 0 also
% compute H+H
daysH2_ = conv(daysHfrac_,daysHfrac_); % note: support = 0:120
% move the (small) support at 0 to 1 since data has no support at 0
daysH2_(1) = 0;
daysH2_(2) = 1-sum(daysH2_(3:end));

% (source: [1])
NobsW = 4039;
daysWfrac = [0.2 0.6 0.8 0.6 0.9 1.1 1.5 2.2 2.9 2.4 3.5 3.6 ...
             3.2 3.5 3.3 3.3 3.0 2.8 2.6 2.5 2.7 2.3 2.8 1.9 ...
             1.9 1.6 1.8 1.5 1.5 1.6 35.9];
daysWfrac = daysWfrac/sum(daysWfrac);
daysW = 1:31; % (note: final quite large fraction is > 30 days)

[phat,phatci] = lognfit(daysW,0.05,[zeros(1,numel(daysW)-1) 1],daysWfrac);
[gammaW,vargammaW] = lognstat(phat(1),phat(2));
gammaW = 1/gammaW; % = 1/31.9!

% similarly impute the empirical distribution to match H+H above
daysW_ = 31:120;
daysWfrac_ = [daysWfrac(1:end-1) lognpdf(daysW_,phat(1),phat(2))];
daysWfrac_ = [0 daysWfrac_/sum(daysWfrac_)];
daysW_ = [0 daysW(1:end-1) daysW_];

% solve as a constrained least square problem with penalty
fun = @(p)norm(l_convrhs(p,daysH2_)-daysWfrac_)^2+0.1*norm(diff(p))^2;
ndof = numel(daysH2_);
p = exp(1/(1/gammaW-2/gammaH)*daysW_)'; p = p/sum(p);
if ill
  opts = optimoptions(@fmincon,'MaxFunctionEvaluations',1e5);
else
  opts = optimoptions(@fmincon,'MaxFunctionEvaluations',1e5, ...
                      'Display','off');
end
p = fmincon(fun,p,[],[],ones(1,ndof),1, ...
            zeros(ndof,1),[1; ones(ndof-1,1)],[], ...
            opts);

% inspect solution
if ill
  figure(3), clf,
  subplot(2,1,1);
  r = l_convrhs(p,daysH2_);
  bar(daysW_,daysWfrac_); hold on
  plot(daysW_(1:ndof),r,'ro-'), xlim([0 61]);
  legend('Data (incl imputed)','Fit to H+W+H');
  xline(31,'k--','Imputation','HandleVisibilit','off','interpreter','latex');
  subplot(2,1,2);
  plot(daysW_(1:ndof),p(1:ndof),'ro-'), xlim([0 61]);
  title('Exit time intensive care');
  legend('Estimate W');
  xlabel('days');
end

% (2.2) now, reconstruct data and do Bayesian inference as before
daysV = daysW_(1:32);
daysVfrac = p(1:32)';
daysVfrac(end) = daysVfrac(end)+1-sum(daysVfrac);

prior = [1 1/14];
lam = expfit(daysV,0.05,[zeros(1,numel(daysV)-1) 1],daysVfrac);
posteriorW = [prior(1)+NobsW 1/(prior(2)+NobsW*lam)];
if ill
  hold on,
  plot(daysV,exppdf(daysV,1./gamstat(posteriorW(1),posteriorW(2))),'k.-', ...
       'DisplayName','Bayesian fit');
  plot(daysV,exppdf(daysV,1./gaminv(0.05/2,posteriorW(1),posteriorW(2))),'k--', ...
       'HandleVisibility','off');
  plot(daysV,exppdf(daysV,1./gaminv(1-0.05/2,posteriorW(1),posteriorW(2))),'k--', ...
       'HandleVisibility','off');
end

% Bayesian prior (incl simple beta-fit) vs current prior
if ill, figure(4), clf, hold on, end
xx = gaminv([5e-2/2 1-5e-2/2],posteriorW(1),posteriorW(2));
xx = linspace(xx(1),xx(2),1e3);
if ill
  plot(1./xx,1./xx.^2.*gampdf(xx,posteriorW(1),posteriorW(2)),'b');
end
% here is the proposed new fit:
m = mean(xx([1 end]));
s = 0.25*diff(xx([1 end]));
xx = linspace(m-8*s,m+8*s);
gammaWPQ = xx([1 end]);
gammaWHyp = [2 2];
if ill
  plot(1./xx,1./xx.^2.*l_scaledbetapdf(xx,gammaWHyp,gammaWPQ),'k--');
end

% current prior
gammaWPQold = [1/15 1/10];
gammaWMU = 1/11.89;
gammaWHypold(1) = 2;
gammaWHypold(2) = gammaWHypold(1)/((gammaWMU-gammaWPQold(1))/ ...
                             diff(gammaWPQold))-gammaWHypold(1);
xx = linspace(gammaWPQold(1),gammaWPQold(2),1e3);
if ill
  plot(1./xx,1./xx.^2.*l_scaledbetapdf(xx,gammaWHypold,gammaWPQold),'r')
  title('Hyperparameter density: intensive care');
  legend('Bayesian $\Gamma$-fit 95\% CI','Proposed beta-fit',['Current ' ...
                      'prior'], 'interpreter','latex');
  xlabel('days');

  % illustration of fit (here: fit to the orginal data)
  figure(5), clf
  % original gamma-fit
% $$$ sample = exprnd(1./gamrnd(posteriorW(1),posteriorW(2),1,1e5))+ ...
% $$$          exprnd(1./gamrnd(posteriorH(1),posteriorH(2),1,1e5))+ ...
% $$$          exprnd(1./gamrnd(posteriorH(1),posteriorH(2),1,1e5));
% $$$ sampledens = histc(sample,[0 daysW(1:end-1) inf]);
% $$$ sampledens = sampledens(1:end-1)/sum(sampledens);

  % relaxed beta-fit
  sample1 = l_scaledbetarnd(gammaWHyp,gammaWPQ,1,1e5);
  sample2 = l_scaledbetarnd(gammaHHyp,gammaHPQ,1,1e5);
  sample3 =  l_scaledbetarnd(gammaHHyp,gammaHPQ,1,1e5);
  sample = exprnd(1./sample1)+exprnd(1./sample2)+exprnd(1./sample3);
  sampledens = histc(sample,[0 daysW(1:end-1) inf]);
  sampledens = sampledens(1:end-1)/sum(sampledens);

  % current prior
  sample1_ = l_scaledbetarnd(gammaWHypold,gammaWPQold,1,1e5);
  sample2_ = l_scaledbetarnd(gammaHHypold,gammaHPQold,1,1e5);
  sample3_ = l_scaledbetarnd(gammaHHypold,gammaHPQold,1,1e5);
  sample_ = exprnd(1./sample1_)+exprnd(1./sample2_)+exprnd(1./sample3_);
  sampledens_ = histc(sample_,[0 daysW(1:end-1) inf]);
  sampledens_ = sampledens_(1:end-1)/sum(sampledens_);

  bar(daysW,[daysWfrac; sampledens; sampledens_]'), hold on
  plot(daysW,lognpdf(daysW,phat(1),phat(2)),'ro-');
  title('Total exit time intensive care');
  legend('Data','Proposed fit','Current prior','lognormal');
  xlabel('days');
  ylim([0 0.1]);
end

end
% ----------------------------------------------------------------------
function rhs = l_convrhs(p,q)
%L_CONVRHS Right hand side for deconvolution problem.

% this is the leftmost part of P+Q when understood as frequency tables
rhs = conv(p,q);
rhs = rhs(1:numel(q));

end
% ----------------------------------------------------------------------
