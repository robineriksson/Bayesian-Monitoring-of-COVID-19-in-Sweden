%%
%SLAM_RUN is the main script the for SLAM procedure.
%   Consider reading SLAM_INIT to se more what is generated when the
%   problem is properly structured.
%
%   main output generated from this script is THETAS which is the
%   approximated posterior.
%
% R.Eriksson 2021/01/21.

region='Uppsala';
nMCMC = 1e2;
smc=0; % SMC-sampler? > 0 means yes, and the number of generations.
verb = true; % debug text
nslab=1; % {1: 1st, 8: 8th, 15: 15th, 22:22nd} monthly slabs
date=[200401 210531];

prefix = postpath;
prefix = [prefix '/SLAM/perRegion/'];
%posterior = 'perRegion/slam201231_22_Uppsala_monthly_8.mat';
%posterior = 'perRegion/slam210307_22_Uppsala_monthly.mat';
%posterior = 'perRegion/slam210314_22_Uppsala_monthly.mat';
%posterior = 'perRegion/slam210321_22_Uppsala_monthly.mat';
%posterior = 'perRegion/slam210329_22_Uppsala_monthly.mat';
%posterior = 'perRegion/slam210404_22_Uppsala_monthly.mat';
%posterior = 'perRegion/slam210409_22_Uppsala_monthly_1.mat';
%posterior = 'perRegion/slam210412_22_Uppsala_monthly_8.mat';
%posterior = 'perRegion/slam210418_22_Stockholm_monthly_15.mat';
%posterior = 'perRegion/slam210418_22_Uppsala_monthly_15.mat';
%posterior = 'perRegion/slam210420_22_Stockholm_monthly_15.mat';
%posterior = 'perRegion/slam210420_22_Uppsala_monthly_15.mat';
%posterior = 'perRegion/slam210420_22_super1_monthly_15.mat';
%posterior = 'perRegion/slam210418_22_Stockholm_monthly_15.mat';
%posterior = 'perRegion/urdme_slam210418_22_Stockholm_monthly_15.mat';
%posterior = 'perRegion/slam210424_22_Stockholm_monthly_22.mat';
%posterior = 'perRegion/slam210425_22_Uppsala_monthly_22.mat';
%posterior = 'perRegion/slam210502_22_Uppsala_monthly_1.mat';
%posterior = 'perRegion/slam210331_Uppsala_monthly_1.mat';
%posterior = 'perRegion/slam210509_Stockholm_monthly_8.mat';
%posterior = 'perRegion/slam210509_Uppsala_monthly_8_100.mat';
%posterior = 'perRegion/slam210517_Uppsala_monthly_15_100.mat';
%posterior = 'perRegion/slam210518_Uppsala_monthly_15_100.mat';
%posterior = 'perRegion/slam210518_Stockholm_monthly_15_100.mat';
%posterior = 'perRegion/slam200601_Stockholm_monthly_8_100.mat';
%posterior = 'perRegion/slam210331_Uppsala_monthly_1_100.mat';
%posterior = 'perRegion/slam200601_Uppsala_monthly_1_100.mat';
%posterior = 'perRegion/slam210526_Uppsala_monthly_22.mat';
%posterior = 'perRegion/slam210526_Stockholm_monthly_22.mat';
%posterior = 'slam210531_Uppsala_monthly_1_URDME1_100.mat';
posterior = 'slam210531_Uppsala_monthly_1_100.mat';

scaleS=0.05; % > 1 if smc.

register='C19';

lsmoothing = false;

useinit   = true; % use old posterior as initalization
useState0 = false; % use an origo file

useCSSS = false; % check this
fix=false; % almost always false


% old flags, that still are needed, but with no effect.
perslab=false;
% initial guess
init = getInit(useinit,[prefix posterior],false);


% use initial state to start at a later date
state0 = getInitState(useState0,210331,region);
  

rng(0)
[thetas, sl, slab, amparam, amfunc, outverb,Ydata] = slam_init(region,nMCMC,...
  verb,nslab,date,init,scaleS,register,useCSSS,perslab,fix,...
  smc,state0);

return;
%% store results
burnin = 1e0;
jump=1e1;
tosave=false;

[rates,rates100,ratesR0] = savePosterior(thetas, sl, slab, ...
  amparam, burnin, jump, ...
  true, useCSSS, tosave, fix, nslab,smc,[],{'full' '100'});
    
return;
%% continue chain
% If the above has been runned, this code can be used to reuse the
% variables already used, and increase the length of the parameter chain.
start = size(thetas,2);
nMCMC2 = 1e5;
stop = start + nMCMC2;

thetas = [thetas, zeros(size(thetas,1), nMCMC2)];
sl = cat(2,sl, zeros(size(sl,1),nMCMC2));

tic;
[thetas, sl, amparam, outverb] = amfunc.metro(nMCMC, thetas, sl, amparam, amfunc, Ydata, start, stop);
tt = toc;

return;
%% plot only R0 trajectories.
figure(1)
jump = 1e1;
burnin= 1e0;
subplot(3,4,1)
plot(thetas(ismember(amparam.ratenames,'R0'),burnin:jump:end)');
title('R0')
subplot(3,4,2)
plot(thetas(ismember(amparam.ratenames,'half_life'),burnin:jump:end)');
title('half life')
subplot(3,4,3)
plot(thetas(ismember(amparam.ratenames,'IFR'),burnin:jump:end)');
title('IFR')
subplot(3,4,4)
plot(thetas(ismember(amparam.ratenames,'IC_HOSP'),burnin:jump:end)');
title('IC-HOSP')
% plot(sum(sl,1));
% title('SL')

subplot(3,4,5)
plot(1./thetas(ismember(amparam.ratenames,'sigma'),burnin:jump:end)');
ylim([2 8])
%yline(1./D.rates.sigma)
title('1/sigma')
subplot(3,4,6)
plot(1./thetas(ismember(amparam.ratenames,'gammaI'),burnin:jump:end)');
%yline(1./D.rates.gammaI)
title('1/gammaI')
subplot(3,4,7)
plot(1./thetas(ismember(amparam.ratenames,'gammaH'),burnin:jump:end)');
%yline(1./D.rates.gammaH)
title('1/gammaH')
subplot(3,4,8)
plot(1./thetas(ismember(amparam.ratenames,'gammaW'),burnin:jump:end)');
%yline(1./D.rates.gammaW)
title('1/gammaW')





subplot(3,4,9)
plot(thetas(ismember(amparam.ratenames,'HOSP'),burnin:jump:end)');
title('HOSP')

subplot(3,4,10)
plot(thetas(ismember(amparam.ratenames,'thetaA_'),burnin:jump:end)');
title('thetaA')
subplot(3,4,11)
plot(thetas(ismember(amparam.ratenames,'thetaE_'),burnin:jump:end)');
title('thetaE')
subplot(3,4,12)
plot(thetas(ismember(amparam.ratenames,'E2I'),burnin:jump:end)');
title('E2I')

return;
%% plot used data - for a single region
figure(2)
semilogy(squeeze(Ydata)')
for i = 1:numel(amparam.slabstop)
  xline(find(amparam.date == amparam.date(amparam.slabstop(i))))
end
x=xticks();
xticklabels(amparam.date(x(1:end-1)+1))
xtickangle(45)