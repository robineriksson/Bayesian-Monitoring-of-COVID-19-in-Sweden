%COVID19ENGER_RUN_POST Runs the EngEr Covid-19 model extending
% COVID19ENGER_RUN by also allowing for posterior noise by sampling the
% posterior further than the mean.

% R. Marin 2022-03-14 (added posterior sampling support)
% S. Engblom 2021-04-17 (Major revision)
% S. Engblom 2020-06-26 (Revision)
% S. Engblom 2020-04-09 (Revision)
% S. Engblom 2020-03-28


rng(200328); % reproducible results

% simulation input
inline = false; % use inline propensities or not

% pick region(s) + date
%reg = 1; % *** region select
%date = '210523';
regionList = regions(false);
Nvoxels = numel(reg);
fprintf('\n*** EngEr C19 in %d region(s): Nreplicas = %d, Psamples = %d ***\n', ...
  Nvoxels,Nreplicas,Psamples);

%---------------------------------------------------------------------------

%% (1) build the model
disp('Building model...');
if inline
  umod = covid19enger_inline;
else
  if divide
    umod = covid19enger2;
  else
    umod = covid19enger;
  end
end
umod.vol = ones(1,Nvoxels);
umod.sd = ones(1,Nvoxels);
Nspecies = size(umod.N,1);
umod.D = sparse(Nspecies*Nvoxels,Nspecies*Nvoxels);
D = sparse(Nvoxels,Nvoxels);
sumD = sum(D,1)';

% posterior rates
path = mfilename('fullpath');
resultdir = postpath();
posterior = cell(1,Nvoxels);
%u0 = []; phi0 = zeros(1,Nvoxels);
for i = 1:Nvoxels
  postfile = ['slam' date '_' regionList{reg(i)} '_monthly_1_100.mat'];
  posterior{i} = strcat([resultdir 'KLAM/perRegion/'],postfile);
end

% weighted posterior sample
weights = umod.private.POP;
weights = weights(reg);
P = posteriorenger(Psamples,posterior,weights);
Nslabs = numel(P.meta.slabstop)-1;
if isinf(Psamples)
  Psamples_ = 1;
else
  Psamples_ = Psamples;
end


% define a common time frame, including dates
DATES = P.meta.date;
DATES = DATES(DATES <= str2double(stopdate));
TSPAN = 1:numel(DATES);
Usamples = zeros(Nspecies*Nvoxels,numel(DATES),Nreplicas,Psamples_);
for p = 1:Psamples_
  fprintf('Post sample: %d',p)
    % URDME fractions and rates
  clear Rates
  Rates.F0 = P.F0ave(:,p);
  Rates.F1 = P.F1ave(:,p);
  Rates.F2 = P.F2ave(:,p);
  Rates.F2d = P.F2dave(:,p);
  Rates.F3 = P.F3ave(:,p);
  Rates.F3d = P.F3dave(:,p);
  Rates.F4 = P.F4ave(:,p);
  Rates.sigma = P.sigma(p);
  Rates.gammaI = P.gammaI(p);
  Rates.gammaA = P.gammaA(p);
  Rates.gammaH = P.gammaH(p);
  Rates.gammaW = P.gammaW(p);

  % handled manually later
  % $$$ Rates.rho = P.rho;
  % $$$ Rates.thetaI = P.thetaI;
  % $$$ Rates.thetaA = P.thetaA;
  % $$$ Rates.thetaE = P.thetaE;
  % $$$ Rates.lambda = 0*P.lambda;

  if inline
    warning('Inline propensities not tested for a while.');
    % inline uses a different sets of rates and they are not fed via the
    % .ldata/.gdata fields
    Rates = inline_rates(Rates);
    GDATA = zeros(numel(umod.private.RateNames)-1,Nslabs);
    for r = 2:numel(umod.private.RateNames) % beta is first, treated separately
      GDATA(r-1,:) = Rates.(umod.private.RateNames{r})';
    end
    GDATA = GDATA(umod.private.ixRate(2:end)-1,:);
  else
    % *** this construct relies on using the same ordering of rates:
    if ~isequal(fieldnames(Rates),umod.private.RateNames(3:end)')
      error('It seems the ordering of rate parameters differ.');
    end

    % convert into Mrates-by-Nslabs matrix:
    GDATA = struct2cell(Rates);
    GDATA = [cat(2,GDATA{1:7})'; ...
      repmat(cell2mat(GDATA(8:end)),1,Nslabs)];
  end

  % handled manually
  rho = P.rho(p);
  thetaI = P.thetaI(p);
  thetaA = P.thetaA(p);
  thetaE = P.thetaE(p);
  lambda = 0*P.lambda(:,p); % note: lambda = 0
  disp('   ...done.');

  %% (2) initial conditions & simulation time interval
  if ~strcmp(solver,'uds')
    disp('Initial data...');
  end


  % initial data from file
  umod.u0 = [u0; zeros(1,Nvoxels)];
  % the 7 states are [I A E phi H W D]; add also state R(ecovered)

  umod.tspan = TSPAN;
  if ~strcmp(solver,'uds')
     disp('   ...done.');
  end


  %% (3) parse and compile, prepare to solve
  if ~strcmp(solver,'uds')
    disp('Parse and compile...');
  end
  umod = urdme(umod,'solve',0,'compile',p==1,'solver',solver, ...
    'gdata',zeros(size(umod.gdata)), ...
    'ldata',zeros(size(umod.ldata,2),Nvoxels));
  EulFwd = strcmp(solver,'uds');
  umod.solve = 1;
  umod.parse = 0;
  umod.compile = 0;

  % solve in split-step fashion
  tspan = umod.tspan;
  U = zeros(Nspecies*Nvoxels,numel(tspan),Nreplicas);
  U(:,1,:) = repmat(umod.u0(:),1,1,Nreplicas);
  PHI = zeros(Nvoxels,numel(tspan),Nreplicas);
  PHI(:,1,:) = repmat(phi0(:),1,1,Nreplicas);
  if ~strcmp(solver,'uds')
    disp('   ...done.');
  end



  %% (4) solve
  if ~strcmp(solver,'uds')
    disp('Solving...');
  end

  for i = 2:numel(tspan)
    slab_ = P.meta.slabs(i-1);
    umod.tspan = tspan(i-1:i);
    if ~strcmp(solver,'uds')
      fprintf('%d%% ',ceil(100*i/numel(tspan)));
    end

    if mod(i,20) == 0, fprintf('\n'); end
    % update local state using present value of phi
    for k = 1:Nreplicas
      umod.u0 = reshape(U(:,i-1,k),Nspecies,Nvoxels);
      if inline
        umod.inline_propensities.K(umod.private.ixK(1)) = beta(i-1);
        umod.inline_propensities.K(umod.private.ixK(2:end)) = GDATA(:,slab_);
      else
        umod.ldata = [beta(i,reg); ...
          PHI(:,i-1,k)'];
        umod.gdata = GDATA(:,slab_);
      end
      umod.seed = randi(intmax);
      if EulFwd
        dU = umod.N*reshape(mexrhs(tspan(i-1),umod.u0,size(umod.N,2), ...
          umod.vol,umod.ldata,umod.gdata,umod.sd,[],[],[]),size(umod.N,2),[]);
        umod.U = U(:,i-1,k)+reshape(dU,[],1);
      else
        umod = urdme(umod);
      end
      U(:,i,k) = umod.U(:,end);
    end

    % update phi
    UI = reshape(U(1:Nspecies:end,i-1,:),Nvoxels,Nreplicas);
    UA = reshape(U(2:Nspecies:end,i-1,:),Nvoxels,Nreplicas);
    UE = reshape(U(3:Nspecies:end,i-1,:),Nvoxels,Nreplicas);
    for k = 1:Nreplicas
      % commuter shedding
      shed = lambda(slab_)*( ...
        ... % incoming/outgoing
        D*(thetaE*UE(:,k)+thetaA*UA(:,k)+thetaI*UI(:,k))- ...
        sumD.*(thetaE*UE(:,k)+thetaA*UA(:,k)+thetaI*UI(:,k)));
      % local shedding
      shed = (shed+thetaE*UE(:,k)+thetaA*UA(:,k)+thetaI*UI(:,k));
      % exponential exact integrator for linear decay rho
      PHI(:,i,k) = PHI(:,i-1,k)*exp(-rho)-shed/rho*expm1(-rho);
    end
    % put phi where it belongs
    U(4:Nspecies:end,i,:) = PHI(:,i,:);
  end
  Usamples(:,:,:,p) = U;
end
umod.tspan = tspan;
umod.U = Usamples;
fprintf('\n');

disp('   ...done.');
