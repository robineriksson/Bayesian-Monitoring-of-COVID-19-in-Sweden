function [x00,beta,J,misfit,beta_mpc,ix_mpc] = mpcloop_beta_S(y,F,H,...
    betaIdx,wbeta,slab,S,steplength,windowlength,x0,beta0,options)
%MPCLOOP_BETA_S solves dynamic beta optimization using an MPC-setup
% Same as MPCLOOP_BETA, except for support for slabs, and covariance matrix
% S which scales the measurements.

% H. Runvik 2021-04-12

if nargin<12
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',1700000,...
        'MaxIterations',6000,'Display','iter-detailed');
end

nstate = size(F{1},1);
ntime = size(y,2);
nstart=[1:steplength:ntime-windowlength ntime-windowlength+1];
nsteps = numel(nstart);
beta=zeros(1,ntime);

betastruct.idx=betaIdx;
beta_mpc = zeros(nsteps,windowlength);
ix_mpc = zeros(nsteps,windowlength);

for k=1:nsteps
    if k==1
        y_loc = y(:,nstart(k):nstart(k)+windowlength-1);
        beta0_loc = beta0(1:windowlength);
%         if max(abs(diff(beta0_loc)))>1e-5
%             ub = beta0_loc.*[1.1*ones(1,floor(windowlength*0.7))...
%                 linspace(2, 3, ceil(windowlength*0.3))];
%             lb = beta0_loc.*[0.9*ones(1,floor(windowlength*0.7))...
%                 linspace(0.5, 0, ceil(windowlength*0.3))];
%         else
            ub = 5*beta0_loc;
            lb = 0*beta0_loc;
%         end
        slab_loc = slab(nstart(k):nstart(k)+windowlength-1);
        betaFunc = @(u)betaCost_S(u,y_loc,F,H,betaIdx,wbeta,slab_loc,S(:,:,nstart(k):nstart(k)+windowlength-1),x0);
        [betatemp,~,exitflag]=fmincon(betaFunc,[beta0_loc],[],[],[],[],...%[x0' beta0_loc]
            lb,...%+nstate
            ub,[],options);%[inf*ones(1,nstate) 0.4*ones(1,windowlength)]
        if exitflag ~=1 && exitflag ~=2
            disp('Convergence problem')
        end
        if min(betatemp-lb)<1e-6 || min(ub-betatemp)<1e-6
            disp('Constraint active')
        end
        %x00=betatemp(1:nstate);
        x00=x0';
    else
        options.InitBarrierParam=1e-4;
        y_loc = y(:,nstart(k):nstart(k)+windowlength-1);
        slab_loc = slab(nstart(k):nstart(k)+windowlength-1);
%         if max(abs(diff(beta0_loc)))>1e-5
%             ub = beta0_loc.*[1.1*ones(1,floor(windowlength*0.7))...
%                 linspace(2, 3, ceil(windowlength*0.3))];
%             lb = beta0_loc.*[0.9*ones(1,floor(windowlength*0.7))...
%                 linspace(0.5, 0, ceil(windowlength*0.3))];
%         else
            ub = 5*beta0_loc;
            lb = 0*beta0_loc;
%         end
        betaFunc = @(u)betaCost_S(u,y_loc,F,H,betaIdx,wbeta,slab_loc,S(:,:,nstart(k):nstart(k)+windowlength-1),x0);
        [betatemp,~,exitflag]=fmincon(betaFunc,beta0_loc,[],[],[1 zeros(1,windowlength-1)],...
            betaPrev,lb,ub,[],...
            options);
        if exitflag ~=1 && exitflag ~=2
            disp('Convergence problem')
        end
        if min(betatemp-lb)<1e-6 || min(ub-betatemp)<1e-6
            disp('Constraint active')
        end
    end
    
    if k==1
        beta(nstart(k):nstart(k)+steplength-1) = ...
            betatemp(1:steplength);%betatemp(nstate+1:nstate+steplength)
        betastruct.vals=betatemp(1:steplength)';%betatemp(nstate+1:nstate+steplength)
        xSim = LPVsim_slabs(F,x0,betastruct,slab_loc);%betatemp(1:nstate)'
        x0 = xSim(:,end);
        betaPrev = betatemp(steplength+1);%betatemp(nstate+steplength+1)
        beta0_loc = [betatemp(steplength+1:end)...%nstate+steplength+1:end
            beta0(windowlength+1:windowlength+steplength)];
        beta_mpc(k,:) = betatemp(1:end);%(nstate+1:end)
        ix_mpc(k,:) = 1:windowlength;
    elseif k<nsteps
        beta(nstart(k):nstart(k+1)-1)...
            = betatemp(1:nstart(k+1)-nstart(k));
        betastruct.vals = betatemp(1:nstart(k+1)-nstart(k))';
        xSim = LPVsim_slabs(F,x0,betastruct,slab_loc);
        x0 = xSim(:,end);
        betaPrev = betatemp(nstart(k+1)-nstart(k)+1);
        beta0_loc = [betatemp(nstart(k+1)-nstart(k)+1:end)...
            beta0(nstart(k)+windowlength:nstart(k+1)+windowlength-1)];
        beta_mpc(k,:)=betatemp;
        ix_mpc(k,:)=nstart(k):nstart(k)+windowlength-1;
    else
        beta(nstart(k):end) = betatemp;
        beta_mpc(k,:)=betatemp;
        ix_mpc(k,:)=nstart(k):nstart(k)+windowlength-1;
    end
    
end
x00=x00';
[J,misfit] = betaCost_S(beta,y,F,H,betaIdx,wbeta,slab,S,x00);