function C19 = getC19syst(rates,m,n)
%GETC19SYST system matrices for Covid-19 model.
%   C19 = GETC19SYST(rates,m,n) assembles, for rate parameter
%   rates(n), (or rates(m,n) for parameters given as a matrix), a
%   struct C19 containing a transition matrix transMat, a
%   deterministic update law dMat, and a list of
%   cumulative/incidence/deterministic/removed states
%   CStates/IStates/DStates/Astates. The output is intended to be used
%   as input to GETC19FILT.
%   
%   C19.transMat(i,j) defines the rate of the transition from the
%   compartment i to compartment j.
%   
%   C19.dMat defines the deterministic linear dynamics, used to
%   describe the evolution of the infectious pressure. It is a matrix,
%   whose rows correspond to the contributions from the states defined
%   by the rows of transMat, to these dynamics, i.e., for a scalar
%   deterministic state phi_k and a state vector x_k, the update law
%   is phi_{k+1} = dMat*x_k.
%   
%   CStates, IStates, DStates, and AStates are arrays of indices which
%   define properties of the corresponding states (given by transMat)
%   in the following way:
%
%   -C19.CStates: states that will get cumulative counterparts in
%   order to, e.g., allow for comparison with measured data.
%
%   -C19.IStates: states that similarly will get incidence
%   counterparts.
%
%   -C19.DStates: deterministic states. The corresponding columns of
%   transMat are all zero, since transitions into deterministic states
%   is not supported.
%
%   -C19.AStates: states that will be removed in the final model. The
%   purpose of including states in transMat that will be removed later
%   is to represent the transitions from other states into these sink
%   states in a convenient way.
%   
%   See also GETC19FILT, C19FILT_KALMAN.

% S. Engblom 2021-03-03 (Revision, optional CStates/IStates)
% S. Engblom 2020-10-22 (Revision, slab)
% H. Runvik 2020-10-01 (Revision, adding P state)
% H. Runvik 2020-09-14

beta = rates.beta(m,n);
sigma = rates.sigma(n);
gammaI = rates.gammaI(n);
gammaA = rates.gammaA(n);
gammaH = rates.gammaH(n);
gammaW = rates.gammaW(n);

% scalar (non-age dependent) fractions this middle step of using m_
% instead of m, allows for legacy rate structs (consider removing if
% future rates only involves multi slab F*)
if size(rates.F0ave,1) == 1
  m_ = 1;
else
  m_ = m;
end

F0 = rates.F0ave(m_,n);
F1 = rates.F1ave(m_,n);
F2 = rates.F2ave(m_,n);
F2d = rates.F2dave(m_,n);
F3 = rates.F3ave(m_,n);
F3d = rates.F3dave(m_,n);
F4 = rates.F4ave(m_,n);

thetaI = rates.thetaI(n);
thetaA = rates.thetaA(n);
thetaE = rates.thetaE(n);

% *** testing static rho
% old: rho = rates.rho(m,n); 
rho = rates.rho(n);
% this scaling uses that thetaI = rho:
erho = -expm1(-rho); % (skip division by rho here)
% (thetaI_ = 1 implicitly understood)
thetaA_ = thetaA/thetaI;
thetaE_ = thetaE/thetaI;

% States I A E phi H W D R
C19.transMat = ...
   [0         0            0    0  gammaI*F2     0         gammaI*F2d gammaI*(1-F2-F2d); ...
    gammaA*F1 0            0    0  0             0         0          gammaA*(1-F1); ...
    sigma*F0  sigma*(1-F0) 0    0  0             0         0          0; ...
    0         0            beta 0  0             0         0          0; ...
    0         0            0    0  0             gammaH*F3 gammaH*F3d gammaH*(1-F3-F3d); ...
    0         0            0    0  gammaW*(1-F4) 0         gammaW*F4  0; ...
    0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 0 0 0];

C19.dMat = [erho thetaA_*erho thetaE_*erho exp(-rho) 0 0 0 0];

% cannot be changed for the C19-model:
C19.DStates = 4; % deterministic compartments
C19.AStates = []; % states that should be removed

% optional cumulative/incidence variables:
if isfield(rates,'CStates')
  C19.CStates = rates.CStates; % compartments which get added cumulative state
else
  C19.CStates = [];
end
if isfield(rates,'IStates')
  C19.IStates = rates.IStates; % compartment which gets added incidence state
else
  C19.IStates = [];
end
