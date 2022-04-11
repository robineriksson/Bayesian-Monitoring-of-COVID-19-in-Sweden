%% Generate the data
rng(0)

nMCMC    = 1e4;
smc      = 0;
verb     = false; % debug text
nslab    = 1; % {1: 1st, 8: 8th, 15: 15th, 22:22nd} monthly slabs
date     = [200401 210531];
scaleS   = 0.04; %0.1 in Uppsala gives acceptance rate ~ 10-20%
register = 'URDME1';

useinit      = true; % use old posterior as initalization
useinitstate = false; % use new origo

useCSSS  = false;
perslab  = false;
fix      = false;





reg = [1]; %1:21;
regionList = regions();



prefix = postpath;
prefix = [prefix '/SLAM/perRegion/'];


burnin = 1e0;
jump=1e0;

%posterior = 'perRegion/old_slam210404_22_Uppsala_monthly.mat';
%posterior = 'perRegion/slam201231_Uppsala_monthly_8.mat';
%posterior = 'perRegion/slam210409_22_Uppsala_monthly_1.mat';
%posterior = 'perRegion/slam210412_22_Uppsala_monthly_8.mat';
%posterior = 'perRegion/old_slam210420_22_Uppsala_monthly_15';
%posterior = 'perRegion/urdme_slam210418_22_Stockholm_monthly_15.mat';
%posterior = 'perRegion/slam210420_22_super1_monthly_15.mat';
%posterior = 'perRegion/slam210424_22_Stockholm_monthly_22.mat';
%posterior = 'perRegion/slam210425_22_Uppsala_monthly_22.mat';
%posterior = 'perRegion/slam210502_22_Uppsala_monthly_1.mat';
%posterior = 'perRegion/slam210331_Uppsala_monthly_1_100.mat';
%posterior = 'perRegion/slam210331_';
%posterior = 'perRegion/slam210331_Östergötland_monthly_1_100.mat';

%posterior = 'perRegion/slam210502_22_Stockholm_monthly_1.mat';
%posterior = 'perRegion/slam210518_';
%posterior = 'perRegion/slam210523_';
%posterior = 'perRegion/slam210509_Stockholm_monthly_8.mat';
%posterior = 'perRegion/slam200601_Stockholm_monthly_8_100.mat';
%posterior = 'perRegion/slam210517_Uppsala_monthly_15_100.mat';
%posterior = 'perRegion/slam210501_Uppsala_monthly_1.mat';
%posterior = 'perRegion/slam210526_Stockholm_monthly_22.mat';
%posterior = 'perRegion/slam210526_Skåne_monthly_1.mat';

tic

round_ = repmat(1:4,1,numel(reg));
reg_ = reshape(repmat(reg,4,1),1,[]);

parpool;
try
  parfor i = 1:numel(round_)
    roundid = round_(i);
    regid = reg_(i);
    region = regionList{regid};
    disp(['starting: ' region ', id: ' num2str(roundid)])

    %posterior = ['slam210531_' region '_monthly_1_100.mat'];
    posterior = ['slam210531_Uppsala_monthly_1_100.mat'];
    if ismember(regid,[1 10 12])
      scaleS = 0.025;
    else
      scaleS = 0.04;
    end

    %  init = getInit(useinit, ...
    %  [prefix posterior region '_monthly_15_100.mat'],false);
    init = getInit(useinit, [prefix posterior],false);
    state0 = getInitState(useinitstate,210331,region);




    [thetas, sl, slab, amparam, amfunc, outverb] = ...
      slam_init(region, nMCMC, verb, nslab, date,...
      init,scaleS,register,useCSSS,perslab,...
      fix,smc,state0);


    disp(['completed: ' region ', id: ' num2str(roundid) ...
      ', early reject rate: ' ...
      num2str(outverb.early_reject_rate) ...
      ', early reject rate 10: ' ...
      num2str(outverb.early_reject_rate_10) ...
      ', early reject rate 10-90: ' ...
      num2str(outverb.early_reject_rate_10_90) ...
      ', early reject rate 90-100: ' ...
      num2str(outverb.early_reject_rate_90_100) ...
      ', accept rate: ' ...
      num2str(outverb.accept_rate)])



    rates = savePosterior(thetas, sl, slab, amparam, burnin, ...
      jump, true, useCSSS, true, fix,nslab,smc,roundid,{'full'});

  end
  stopped = 0;
catch
  disp(['*** CATCH ***']);


  p = gcp; delete(p);
  stopped = true;
end
if ~stopped
  p = gcp; delete(p);
end
disp('*** DONE ***');
toc
return;

%% fast eval of multi-chain rates
% date(end) = 210331;
% nslab=1;
burnin=1e3;%1.5e4;
for i = reg
  region=regionList{i}
  try
    for roundid = 1:4
      file = [prefix '/slam' num2str(date(end)) '_' ...
              region '_monthly_' num2str(nslab) '_run' ...
              num2str(roundid)];
      if contains(register,'URDME')
        file = [file '_' register];
      end
      load(file,'rates','amparam');
      [xrates, names] = struct2mat(rates);
      xrates = xrates(ismember(names,amparam.ratenames),:);
      ratenames = names(ismember(names,amparam.ratenames));
      if roundid == 1
        X = xrates(:,burnin:end);
      else
        X = cat(3,X,xrates(:,burnin:end));
      end
    end

    subset = {'sigma' 'gammaI' 'gammaH' 'gammaW' 'thetaE_'};
    for nid = 1:numel(subset)
      ix = ismember(ratenames,subset(nid));
      subplot(numel(subset),1,nid)
      if strcmp(subset{nid}, 'thetaE_')
          plot(squeeze(X(ix,:,:)))
      else
          plot(1./squeeze(X(ix,:,:)))
      end
      title(subset{nid})
    end

    % Gelman-Rubin statistic
    % chain mean
    mx = squeeze(mean(X,2));
    L = size(X,2);
    % grand mean
    mmx = mean(mx,2);
    J = size(X,3);
    % between chain variance
    B = L/(J-1) * sum( (mx - mmx).^2,2);

    % within chain variance
    s = squeeze(var(X,[],2));
    W = mean(s,2);

    % GR-statistic
    R = (W*(L-1)/L + B/L) ./ W;
    R_mean = mean(R)
  catch
    warning(['missing region: ' region])
  end
end
return;



%% Combine all the files, and produce a posterior per region.
% remove-burnin
date_ = num2str(date(end));
ending = num2str(nslab);
burnin = 1e3;
saveit=false; % false if testing

for i = reg
  region = regionList{i};
  try
    post0 = [prefix '/slam' date_ '_' region '_monthly_' ending '_run'];
    file = {[post0 '1'] [post0 '2'] [post0 '3'] [post0 '4']};
    if contains(register,'URDME')
      for i = 1:numel(file)
        file{i} = [file{i} '_' register];
      end
    end


    for i = 1:numel(file)
      try
        load(file{i},'rates','sl_burn');
        [m,n] = struct2mat(rates);
        if i == 1
          mat = m(:,burnin+1:end);
          sl = sl_burn(:,burnin+1:end);
        else
          mat = cat(2,mat,m(:,burnin+1:end));
          sl = cat(2,sl,sl_burn(:,burnin+1:end));
        end
      catch
        warning(['did not find: ' file{i}])
      end
    end

    load(file{1},'amparam','slabs');
    [rates,rates100,ratesR0] = savePosterior(mat, sl, slabs, amparam, 1e0, ...
      1e0, true, useCSSS, saveit, fix,nslab,smc,[],{'full' '100'});
  catch
    warning(['missing region: ' region])
  end
end

% *** old version
% %% construct bagged posterior
% %date = '210429';
% date_ = num2str(date(end));
% ending = '1';
% region = regionList{reg}
% post0 = [prefix 'perRegion/slam' date_ '_22_' region '_monthly_' ending '_run'];
% file = {[post0 '1.mat'] [post0 '2.mat'] [post0 '3.mat'] [post0 '4.mat']};
% rates = posteriorenger(1e3,file,[1 1 1 1]);
% load([post0 '1.mat'],'amparam','slab');
% rates.meta.hash_all = rates.meta.hash;
% rates.meta.hash = fsetop('check', rates.meta.hash'); % hash of hashes
% outfile = [prefix 'perRegion/slam' date_ '_22_' region '_monthly_' ending '.mat'];
% %save(outfile,'rates','amparam','slab');
% return;