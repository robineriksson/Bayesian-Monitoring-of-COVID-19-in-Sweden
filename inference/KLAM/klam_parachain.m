% KLAM_PARACHAIN generates multiple (parallel) Markov Chains using
% Kalman marginal Adaptive Metropolis, see AM.m for the Metropolis
% implementation. Per parallel core, the inference object is
% initialized with KLAM_INIT.m, and then run for NMCMC samples using
% AM.m
%
% After the sampling is done, the parallel chains are evaluated by
% Gelman-Rubin convergence score if GELMRUB = [true]. Each chain is
% saved as a temporary file, and can be recovered after the fact even
% if SAVETOFILE = [false]. But if SAVETOFILE = true, after the
% Gelman-Rubin score is calculated, a posterior file is compiled and
% saved under /inference/results/KLAM/perRegion/.
%
% Adjustable parameters are
% NMCMC = 5e4: length of each individual chain
% SAVETOFILE = false: if the posterior file is to be saved to disk
% REG = 2: the region to compute for, see regions() for the numbering.
% GELMANRUB = true: computes the Gelman-Rubin convergence score,
% betweem 0.9-1.2 is good.
% BURNIN = 1e3: the number of samples to discard in the start of each chain.
% RUNS = 4: number of parallel chains, should be set to the number of
% physical core available for optimal efficiency.
% REGISTER = 'C19': the data  source to use for inference, see LOADDATA(register).

% R. Marin 22-04-17

%% Generate the data
rng(0)

if ~exist('nMCMC','var')
    nMCMC=5e4;
end

if ~exist('savetofile','var')
    savetofile=false;
end

if ~exist('reg','var')
    reg=2;
else
    if max(reg) > 21
        error(['Only supported to run 21 regions,'...
               'National posterior is sampled from a basket of region']);
    end
end

if ~exist('gelmanRub','var')
    gelmanRub=true;
end

if ~exist('burnin','var')
    burnin=1e3;
end

if ~exist('runs','var')
    runs=4;
end

if ~exist('register','var')
    register = 'C19';
end

if ~exist('verb','var')
    verb = false;
end

if ~exist('waste','var')
    waste = false;
end

if ~exist('scaleS','var')
    scaleS = 0.04; % adjust to find good acceptance rate
                   %scaleS = 0.04*8 for Stockholm
end

if ~exist('par','var')
    par=true; % if chains in parallel or not
end

if ~exist('date','var')
    date = [200401 210531]; % what date to look at
end

if ~exist('postdate','var')
    postdate=210531; % what previous posterior to look at
                     % postdate = NaN; if no inital guess.
end


nslab    = 1; % {1: 1st, 8: 8th, 15: 15th, 22:22nd} monthly slabs

% use old posterior as initalization
useinit  = ~isnan(postdate);

% do not change these
useinitstate = false; % use new origo
useCSSS  = false;
perslab  = false;
fix      = false;


regionList = regions(false);



prefix = postpath;
if ispc
    prefix = [prefix 'KLAM\perRegion\'];
else
    prefix = [prefix 'KLAM/perRegion/'];
end


burnin = 1e0;
jump=1e0;




tic
rates_ = l_samplechains_par(par, scaleS, runs, reg, regionList,...
                            useinit, useinitstate,  prefix, nMCMC, ...
                            verb, jump, burnin, nslab, date, ...
                            register, useCSSS, perslab, fix,waste, ...
                            postdate);
toctime = toc;
if verb
    toctime
end

%% fast eval of multi-chain rates

if gelmanRub
    for i = reg
        region=regionList{i};
        if verb
            disp(['region: ' region]);
        end
        try
            file = cell(1,runs);
            for roundid = 1:runs
                if waste
                    file{roundid} = [prefix 'slam' num2str(date(end)) '_' ...
                                     region '_monthly_' num2str(nslab) '_waste_run' ...
                                     num2str(roundid)];
                else
                    file{roundid} = [prefix 'slam' num2str(date(end)) '_' ...
                                     region '_monthly_' num2str(nslab) '_run' ...
                                     num2str(roundid)];
                end
                if contains(register,'URDME')
                    file{roundid} = [filefile{roundid} '_' register];
                end
            end

            for roundid = 1:runs

                load(file{roundid},'rates','amparam');
                [xrates, names] = struct2mat(rates);
                xrates = xrates(ismember(names,amparam.ratenames),:);
                ratenames = names(ismember(names,amparam.ratenames));
                if roundid == 1
                    X = xrates(:,burnin:end);
                else
                    X = cat(3,X,xrates(:,burnin:end));
                end
            end

            %Uncomment below  if you want to see some example trajectories
            %
            subset = {'sigma' 'gammaI' 'gammaH' 'gammaW'};
            for nid = 1:numel(subset)
                 ix = ismember(ratenames,subset(nid));
                 subplot(numel(subset),1,nid)
                 plot(1./squeeze(X(ix,:,:)))
                 title(ratenames{nid})
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
            R_mean = mean(R);
            if verb
                disp(['GR: ' num2str(R_mean)]);
            end
        catch
            error(['GR missing region: ' region])
        end
    end
end


%% Combine all the files, and produce a posterior per region.
% remove-burnin
date_ = num2str(date(end));
ending = num2str(nslab);
% types = {'full' '100'};
types = {'full'};
% types = {'100'};

for i = reg
    region = regionList{i};
    try
        file = cell(1,runs);
        for roundid = 1:runs
            if waste
                file{roundid} = [prefix 'slam' num2str(date(end)) '_' ...
                                 region '_monthly_' num2str(nslab) '_waste_run' ...
                                 num2str(roundid)];
            else
                file{roundid} = [prefix 'slam' num2str(date(end)) '_' ...
                                 region '_monthly_' num2str(nslab) '_run' ...
                                 num2str(roundid)];
            end
            if contains(register,'URDME')
                file{roundid} = [filefile{roundid} '_' register];
            end
        end

        for i = 1:numel(file)
            try
                load(file{i},'rates','sl_burn','slabs');
                [m,n] = struct2mat(rates);
                if i == 1
                    mat = m(:,burnin+1:end);
                    sl = sl_burn(:,burnin+1:end);
                else
                    mat = cat(2,mat,m(:,burnin+1:end));
                    sl = cat(2,sl,sl_burn(:,burnin+1:end));
                end
            catch
                error(['did not find: ' file{i}])
            end
        end

        [rates,rates100,ratesR0] = savePosterior(mat, sl, slabs, amparam, 1e0, ...
                                                 1e0, true, useCSSS, ...
                                                 savetofile, ...
                                                 fix,nslab,0,[],types,false,...
                                                 verb,waste);
    catch
        error(['missing region: ' region ', file: ' file{i}])
    end
end

return;

function c_rates = l_samplechains_par(par, scaleS, runs, reg, regionList,...
                                      useinit, useinitstate,  prefix, nMCMC, ...
                                      verb, jump, burnin, nslab, date, ...
                                      register, useCSSS, perslab, fix,waste,...
                                      postdate)

% L_SAMPLECHAINS_PAR generates RUNS number of chain that can be
% combined. If PAR = TRUE, then the chains are sampled in
% parallel.

    round_ = repmat(1:runs,1,numel(reg));
    reg_ = reshape(repmat(reg,runs,1),1,[]);

    c_rates = cell(numel(round_),1);

    if par
        parpool;
        try
            parfor i = 1:numel(round_)
                c_rates{i} = l_samplechain(i, scaleS, round_, reg_, regionList,...
                                           useinit, useinitstate,  prefix, nMCMC, ...
                                           verb, jump, burnin, nslab, date, ...
                                           register, useCSSS, perslab, fix,waste,...
                                           postdate);
            end
            stopped = 0;
        catch
            p = gcp; delete(p);
            disp('*CATCHED*');
            stopped = true;
        end
        if ~stopped
            p = gcp; delete(p);
        end
        if ~verb
            warning('off','klam:slabs');
        end
    else
        for i = 1:numel(round_)
             c_rates{i} = l_samplechain(i, scaleS, round_, reg_, regionList,...
                                        useinit, useinitstate,  prefix, nMCMC, ...
                                        verb, jump, burnin, nslab, date, ...
                                        register, useCSSS, perslab, fix,waste,...
                                        postdate);
        end
    end
    c_rates = reshape(c_rates,[],numel(unique(reg_)));
end


function rates = l_samplechain(i, scaleS, round_, reg_, regionList,...
                               useinit, useinitstate,  prefix, nMCMC, ...
                               verb, jump, burnin, nslab, date, ...
                               register, useCSSS, perslab, fix,waste,...
                               postdate)

    % L_SAMPLECHAIN generates one sample chain, and is designed to be
    % using inside of a loop, possibly parallel.

    roundid = round_(i);
    regid = reg_(i);
    region = regionList{regid};

    if verb
        disp(['starting: ' region ', id: ' num2str(roundid)])
    end

    ending = num2str(nslab);

    posterior = ['slam' num2str(postdate) '_' region '_monthly_' ending];

    if waste
      posterior = [posterior '_waste'];
    end

    if contains(register,'URDME')
      posterior = [posterior '_' register];
    end

    posterior = [posterior '_100.mat'];


    if ismember(regid,[1 10 12])
        scaleS_ = scaleS*0.25;
    else
        scaleS_ = scaleS*1.00;
    end

    init = getInit(useinit, [prefix posterior],false);
    state0 = getInitState(useinitstate,210331,region);

    if ~verb
        warning('off','klam:slabs');
    end

    [thetas, sl, slab, amparam, amfunc, outverb] = ...
        klam_init(region, nMCMC, verb, nslab, date,...
                  init,scaleS_,register,useCSSS,perslab,...
                  fix,0,state0, waste);

    if verb
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
    end
    rates = savePosterior(thetas, sl, slab, amparam, burnin, ...
                          jump, true, useCSSS, true, fix,nslab,...
                          0,roundid,{'full'},false,verb,waste);

end
