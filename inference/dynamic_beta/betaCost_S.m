function [J, L] = betaCost_S(beta,y,F,H,betaIdx,wbeta,slab,S,x0)
%BETACOST_S cost-function for state space model with time-varying parameter
%beta.
%   [J, L] = BETACOST_S(beta,y,F,H,betaIdx,wbeta,slab,S,x0) calculates cost
%   J and likelihood L (which is scaled and shifted, see below) for data y
%   for the state space model 
%   x_{k+1} = F_k x_k, x_0 = x0,
%   y_k=H x_k
%   where the elements of F_k are constant and given by F, except for the
%   element F(betaIdx), which for each time step k is defined by the
%   corresponding element in the vector beta. The covariance matrix S
%   scales the measurements. J relates to L as
%   J = wbeta*sum(diff(beta).^2) - L,
%   i.e. wbeta penalizes acts as a regularization that penalizes fast
%   variations in beta. Note that the cost is rescaled with a factor 1e-3
%   and that the fixed offsets of the likelihood are not calculated for
%   performance reasons. See MPCLOOP_BETA_S for the use of the function in
%   an optimization formulation.

% H. Runvik 2021-01-12

nstate=size(F{1},1);
ntime=size(y,2);
if nargin<9
    x0=beta(1:nstate)';
    beta=beta(nstate+1:end);
end
L=0;
x=x0;

slabloop = 1+[0 find(diff(slab) ~= 0) numel(slab)];
slabloop = [slabloop(slabloop<=ntime) ntime+1];
nslab = numel(slabloop)-1;
for m = 1:nslab
  Fslab = F{m+slab(1)-1};
  for n=slabloop(m):slabloop(m+1)-1
    
    % exceptional: check for NaN in measurement
    yNum = y(:,n);

    Hmod = H;
    if any(isnan(yNum))
      ixn = find(~isnan(yNum));
      Z = sparse(ixn,ixn,1,numel(yNum),numel(yNum));
      % remove them:
      yNum = Z*yNum;
      Hmod = Z*H;
    end
    
    Fslab(betaIdx)=beta(n);
    x=Fslab*x;
    yy = yNum - Hmod*x;%y(:,n)-H*x;
    L = L -0.5*(yy'*(S(:,:,n)\yy));%+size(y,1)*log(2*pi)+log(det(S(:,:,n))));
  end
end

J = 1e-4*(wbeta*sum(diff(beta).^2) - L);
