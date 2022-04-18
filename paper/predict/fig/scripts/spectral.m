% SPECTRALANALYSISF Analyze the operator F "spectrally", ~ODE solution
% without simulation noise. With this formulation, we can estimate the
% CFR from I, H, and W. The output is in CrI and the comparative total
% dead sources as presented from Socialstyrelsen.

% R. Marin 2022-03-22 (stream lined for the results we want in paper)
% S. Engblom 2021-03-11


%% load posterior
% load posterior
Nsample = 1e3;
path = postpath;
prefix = 'KLAM/perRegion/';
% "super-posterior"
date = '210531';
ending = '_1_100';
regionList = regions(false);
reg = 1:21;
regname = regionList(reg);
files = strcat('slam',date,'_',regname,'_monthly', ending, '.mat');
load Ncounties;


filename = strcat(path,prefix,files);

Weights = sum(N(:,reg));
rates = posteriorenger(Nsample,filename,Weights,true);




%%
rng(222);


k = 1:Nsample;


% The posterior estimate needs to be weighted as the posterior is 'slab
% dependent'.
scaleRepo = 'C19'; % Alt. see FHM.
Data_scale = loadData(scaleRepo);
Data_scale = polishData(Data_scale,'D','Dinc',1);
%Data_scale = smoothData(Data_scale,{'D' 'H' 'W'},{'Dinc' [] []});


% First find which slabs should be exluded
stopdate=210323;


if false
  zeroweight = rates.meta.date(rates.meta.slabstop) > stopdate;
  fullweight = find(1-zeroweight);
  zeroweight = find(zeroweight);

  % however, stopdate does probably not end on a slabstop date, find out how
  % many days on that final slab which should be counted and use is as a
  % weight instead of 0/1.
  finalslab = rates.meta.date(rates.meta.slabstop(fullweight(end)):rates.meta.slabstop(zeroweight(1)));
  finalslabw = numel(find(finalslab <= stopdate)) / numel(finalslab);

  % 0 or 1 (or other) weight multiple.
  weightmultiple = ones(1,numel(rates.meta.slabstop)-1);
  weightmultiple(zeroweight(2:end)-1) = 0;
  weightmultiple(zeroweight(1)-1) = finalslabw;

  % Get the data that the posterior used.
  tspan_post = find(Data_scale.date==rates.meta.date(1)):find(Data_scale.date==rates.meta.date(end));
  dates = Data_scale.date(tspan_post);

  % Count the incidence death per slab
  Dinc = sum(Data_scale.Dinc(tspan_post,:), 2);
  Dinc_cum = zeros(1,numel(rates.meta.slabstop)-1);
  for i = 2:numel(rates.meta.slabstop)
    Dinc_cum(i-1) = sum(Dinc(rates.meta.slabstop(i-1):rates.meta.slabstop(i)));
  end

  % Give zero weight to excluded slabs, and partial weight to partial slabs.
  Dinc_cum = Dinc_cum .* weightmultiple;
  weights = Dinc_cum/sum(Dinc_cum); % normalize
else
  weights = ones(1,numel(rates.meta.date(rates.meta.slabstop) > stopdate)-1);
  weights = weights./sum(weights);
end
%% Case fatality risk - risk of dying when in specific compartment

% Find the median, 68, and 95% CrI
quant = [0.05/2 0.32/2 0.5 1-0.32/2 1-0.05/2];

% I --> D *eventually*
ID = 100*(rates.F2dave+rates.F2ave.*(rates.F3dave+rates.F3ave.*rates.F4ave))./(1-rates.F3ave.*(1-rates.F4ave));
ID = weights * ID;
ID_q = quantile(ID, quant);
CFR_I = l_formatCI(round(ID_q,2,'significant'))

% H --> D *eventually*
HD = 100*(rates.F3dave + rates.F3ave.*rates.F4ave)./(1-rates.F3ave.*(1-rates.F4ave));
HD = weights * HD;
HD_q = quantile(HD, quant);
CFR_H = l_formatCI(round(HD_q,2,'significant'))

% W --> D *eventually*
WD = 100*(rates.F4ave+(1-rates.F4ave).*rates.F3dave) ./ ( 1-rates.F3ave.*(1-rates.F4ave));
WD = weights * WD;
WD_q = quantile(WD, quant);
CFR_W = l_formatCI(round(WD_q,2,'significant'))




%% Distribution of the dead (from what compartment)
% I --> [D_I, D_H, D_W] *eventually*
D_pred_WHI0 = [rates.F2ave.*rates.F3ave.*rates.F4ave./(1-rates.F3ave.*(1-rates.F4ave)); ...
  rates.F2ave.*rates.F3dave./(1-rates.F3ave.*(1-rates.F4ave)); ...
  rates.F2dave];

D_pred_IHW = flipud(cumsum(D_pred_WHI0,1));
D_pred_IHW = 100*D_pred_IHW./sum(D_pred_IHW,1);
D_pred_IHW_q =quantile(D_pred_IHW,quant,2);

%%
% according to
%
% https://www.socialstyrelsen.se/globalassets/1-globalt/covid-19-statistik/statistik-om-slutenvard-av-patienter-med-covid-19/statistik-covid19-vardforlopp.xlsx
% with date 2021-03-23, we have
D_SS_IHW = [NaN 3905 890]';
% https://www.socialstyrelsen.se/globalassets/sharepoint-dokument/dokument-webb/statistik/statistik-covid19-vardforlopp.xlsx
% with date 2022-03-10, we have
% D_SS_IHW = [NaN 6800 2063]';

% SS claims in
% https://www.socialstyrelsen.se/globalassets/sharepoint-dokument/dokument-webb/statistik/statistik-covid19-avlidna.xlsx
% The following total deaths:
Data_SSmort = loadData('SS_swe_mort');

D_SS_tot = nansum(Data_SSmort.Dinc(Data_SSmort.date <= stopdate));
D_SS_IHW(1) = D_SS_tot-nansum(D_SS_IHW); % imputed value

% Now scale the SS values such that we end up with the same total as the
% filtered data. Assuming uniform scaling accross all dead compartments.
D_scale_tot = sum(Data_scale.D(Data_scale.date == stopdate,:));
D_scale = D_scale_tot/D_SS_tot;
D_SS_IHW = D_SS_IHW*D_scale;
disp('number of dead')
round(D_SS_IHW,0)
disp('... of which the total is')
sum(round(D_SS_IHW,0))
disp(['... and scaled uniformly by: ' num2str(D_scale)])


D_SS_IHW_norm = 100*D_SS_IHW/sum(D_SS_IHW);

return;
%%
plot(D_pred_IHW(1,:),D_pred_IHW(2,:),'.')
hold on
plot(D_pred_IHW(1,:),D_pred_IHW(3,:),'.')
plot(D_pred_IHW(2,:),D_pred_IHW(3,:),'.')

axis equal

hold off

return;
%% Validate the rate results with the F^\inty * x0 result.
%   C19 = getC19syst(rates,1,k); % *** bug: both m and n required?
%   KF = getC19filt(C19);
%
%   % I A E phi H W D R
%   % 1 2 3 4   5 6 7 8
%   F = KF.F;
%   F_all = cat(3,F_all,F);
%   % D and R are the only stationary states:
%   % $$$ [V,lam] = eig(F,'vector');
% % $$$ eq = find(abs(lam-1) < eps(1e4));
% % $$$ V(:,eq);
% % when R0 = 1 there are 3 eigenvalues = 1, but only 2 are linearly
% % independent eigenvectors
% % for knowledge of setpoints used
% %abspath = mfilename('fullpath');
% load('c19prior.mat','c19prior')
%
%
% % perform the analysis on various small modifications F_ of F
% F_ = F;
%
% % HWfac: W/H at approximate steady-state
% x0 = zeros(size(F_,2),1);

% H --> D eventually
% F_ = F;
% F_(3,:) = 0;
% F_(:,3) = 0;
% x0(:) = 0;
% x0(5) = 100; % 100% @ H at time = 0
% %F_([5 7],6) = 0; % remove W --> {H,D} ("no 2nd chance at death in H")
% xinf = (F_^256*x0);
% HOSP_MORT = [xinf(end-1) rates.F3dave(k)*100 c19prior.SIR_MORT*c19prior.HOSP_MORT*100];
% HOSP_MORT_all = cat(1,HOSP_MORT,HOSP_MORT_all);
% HD(1) = xinf(end-1);

% W --> D eventually
% F_ = F;
% F_(3,:) = 0;
% F_(:,3) = 0;
% x0(:) = 0;
% x0(6) = 100; % 100% @ H at time = 0
% %F_([5 7],6) = 0; % remove W --> {H,D} ("no 2nd chance at death in H")
% xinf = (F_^256*x0);
% HOSP_MORT = [xinf(end-1) rates.F3dave(k)*100 c19prior.SIR_MORT*c19prior.HOSP_MORT*100];
% HOSP_MORT_all = cat(1,HOSP_MORT,HOSP_MORT_all);
%WD(1) = xinf(end-1)

%
% % IFR: I --> D *eventually*
% F_ = F;
% F_(3,:) = 0;
% F_(:,3) = 0;
% x0(:) = 0;
% x0(1) = 100; % 100% @ I at time = 0
% xinf = (F_^256*x0);
% IFR = [xinf(end-1); ...
%        100*(rates.F2dave(k)+rates.F2ave(k)*(rates.F3dave(k)+rates.F3ave(k)*rates.F4ave(k))/(1-rates.F3ave(k)*(1-rates.F4ave(k))))]
%
% IFR_all = cat(1,IFR_all,IFR);


% % IFR: I --> [D_I,D_H,D_W] *eventually*
% F_ = F;
% F_(3,:) = 0;
% F_(:,3) = 0;
% % I A E phi H W D_I R D_H D_W
% % 1 2 3 4   5 6 7   8 9   10
% F_ = [[F_ zeros(size(F_,1),2)]; zeros(2,2+size(F_,2))];
% F_(9,[5 9]) = [F_(7,5) 1]; F_(7,5) = 0;  % re-route H --> D_H
% F_(10,[6 10]) = [F_(7,6) 1]; F_(7,6) = 0; % ditto W --> D_W
% x0(:) = 0;
% x0(1) = 100; % 100% @ I at time = 0
% %x0(5) = 100; % 100% @ H at time = 0
% x0 = [x0; 0; 0];
% xinf = (F_^256*x0);
%
% % percentage of deaths at [I H W]
% D_IHW = xinf([7 9 10])/sum(xinf([7 9 10]))*100;
%
% % prediction of the above directly from rates
%
% D_pred_IHW = [rates.F2dave(k) ...
%               rates.F2ave(k)*rates.F3dave(k)/(1-rates.F3ave(k)*(1-rates.F4ave(k))) ...
%               rates.F2ave(k)*rates.F3ave(k)*rates.F4ave(k)/(1-rates.F3ave(k)*(1-rates.F4ave(k)))]';
% D_pred_IHW = 100*D_pred_IHW/sum(D_pred_IHW)
% D_pred_IHW_all = cat(1,D_pred_IHW_all,D_pred_IHW');

function [CI] = l_formatCI(X)
%L_FORMATCI converts numeric array X into a string \CIedge{...} format ready
%for LaTeX use.

CI = cell(size(X,1),1);
for k = 1:size(X,1)
  CI{k} = ['\CIedge' '{' num2str(X(k,1)) '}'...
    '{' num2str(X(k,2)) '}'...
    '{' num2str(X(k,3)) '}'...
    '{' num2str(X(k,4)) '}'...
    '{' num2str(X(k,5)) '}'];

end
end