function [L, state, Lsum] = logKL(rates,state0,G,D,obsrates,Q,Ydata,...
  slab,idTime,exception)
% [LSUM, STATE]= LOGKL(RATES,sate0, G,D,OBSRATES,YDATA,SLAB,NTIME)
%   give the marginal likelihood for RATES per SLAB using the Kalman
%   filter in the data in YDATA. Supplying the struct EXCEPTION allows
%   for handeling of the early rejections. L is the likelihood as
%   given by Kalman, per day, and Lsum is the sum per SLAB, and state
%   the final state of system.
%
%   *input*
%   RATES     - model parameters
%   STATE0    - struct with fields Z0 and COVZ0, used as inputs in C19FILT
%   G         - Observation matrix
%   D         - Network matrix (often zeros)
%   OBSRATES  - Kalman filter observation input
%   Q         - Kalman filter uncertainty parameter
%   YDATA     - Observational data
%   SLAB      - Slab indicator
%   IDTIME    - what time steps to consider, e.g., all timesteps in slab #1
%   EXCEPTION - struct {'LB' 'UB' 'SDFAC' 'ABSMAGN'}
%
%   see also C19FILT for futher definition.

% R. Eriksson 210522

if nargin < 10
  exception = [];
end
state = {};
slab_ = slab(idTime);
Ydata_ = Ydata(:,:,idTime);
[state.Z,state.covZ,~,L] = C19filt(state0.Y0,state0.Y0Cov,G,rates,D, ...
  obsrates,Ydata_,slab_,numel(idTime), [], Q, exception);
if nargout > 2
  slabidx = [1,1+find(diff(slab_)),numel(slab_)];
  Lsum=zeros(numel(slabidx)-1,numel(rates.sigma));
  for i = 2:numel(slabidx)
    Lsum(i-1,:) = sum(L(slabidx(i-1):slabidx(i),:),1);
  end
else
  Lsum = [];
end
end