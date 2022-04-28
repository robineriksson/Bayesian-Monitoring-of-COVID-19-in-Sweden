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

verb     = false; % verb text
nslab    = 1; % {1: 1st, 8: 8th, 15: 15th, 22:22nd} monthly slabs
date     = [200401 210531];
scaleS   = 0.04; %0.1 in Uppsala gives acceptance rate ~ 10-20%
useinit      = true; % use old posterior as initalization



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

round_ = repmat(1:runs,1,numel(reg));
reg_ = reshape(repmat(reg,runs,1),1,[]);

parpool;
try
    parfor i = 1:numel(round_)

            %for i = 1:numel(round_)
        roundid = round_(i);
        regid = reg_(i);
        region = regionList{regid};
        if verb
            disp(['starting: ' region ', id: ' num2str(roundid)])
        end

        posterior = ['slam210531_' region '_monthly_1_100.mat'];
        if ismember(regid,[1 10 12])
            scaleS = 0.025;
        else
            scaleS = 0.04;
        end

        init = getInit(useinit, [prefix posterior],false);
        state0 = getInitState(useinitstate,210331,region);


        if ~verb
            warning('off','klam:slabs');
        end
        [thetas, sl, slab, amparam, amfunc, outverb] = ...
            klam_init(region, nMCMC, verb, nslab, date,...
                      init,scaleS,register,useCSSS,perslab,...
                      fix,0,state0);

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
                              jump, true, useCSSS, true, fix,nslab,0,roundid,{'full'},false);

    end
    stopped = 0;
catch
    if verb
        disp(['*** CATCH ***']);
    end


    p = gcp; delete(p);
    stopped = true;
end
if ~stopped
    p = gcp; delete(p);
end
if ~verb
    warning('off','klam:slabs');
end
if verb
    disp('*** DONE ***');
end
toctime = toc;
if verb
    toctime
end

%% fast eval of multi-chain rates
if gelmanRub
    % date(end) = 210331;
    % nslab=1;
    for i = reg
        region=regionList{i};
        if verb
            disp(['region: ' region]);
        end
        try
            for roundid = 1:runs
                file = [prefix 'slam' num2str(date(end)) '_' ...
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

            % Uncomment below  if you want to see some example trajectories
            %
            % subset = {'sigma' 'gammaI' 'gammaH' 'gammaW'};
            % for nid = 1:numel(subset)
            %     ix = ismember(ratenames,subset(nid));
            %     subplot(numel(subset),1,nid)
            %     plot(1./squeeze(X(ix,:,:)))
            %     title(ratenames{nid})
            % end

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
                disp(['GR: ' R_mean]);
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


for i = reg
    region = regionList{i};
    try
        post0 = [prefix 'slam' date_ '_' region '_monthly_' ending '_run'];
        file = strcat(post0, strsplit(num2str(1:runs)));
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
                error(['did not find: ' file{i}])
            end
        end

        load(file{1},'amparam','slabs');
        [rates,rates100,ratesR0] = savePosterior(mat, sl, slabs, amparam, 1e0, ...
                                                 1e0, true, useCSSS, savetofile, ...
                                                 fix,nslab,0,[],{'full' '100'},false);
    catch
        error(['missing region: ' region])
    end
end


return;