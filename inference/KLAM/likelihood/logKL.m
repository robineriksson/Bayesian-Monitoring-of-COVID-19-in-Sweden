function [L, state, Lsum] = logKL(rates,state0,G,D,obsrates,Q,Ydata,...
  slab,idTime,exception)
%[LSUM, STATE]= LOGKL(RATES,sate0, G,D,OBSRATES,YDATA,SLAB,NTIME) give the marginal
%   likelihood for RATES per SLAB using the Kalman filter in the data in
%   YDATA. Supplying the struct EXCEPTION allows for handeling of the early
%   rejections.
%
%   rates     - model parameters 
%   state0    - struct with fields Z0 and COVZ0, used as inputs in C19FILT
%   idTime    - what time steps to consider, e.g., all timesteps in slab #1
%   exception - struct {'LB' 'UB' 'SDFAC' 'ABSMAGN'}
%   LSUM has size (unique(slab),numel(rates)), where numel(rates) :=
%   numel(rates.sigma)
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