function [H,R0,rdiag,R0_p2, rdiag_p2] = getC19obs(obsrates)
%GETC19OBS Observation model for Covid-19 model.
%   [H,R0,rdiag] = GETC19OBS(obsrates) produces an observation matrix
%   H and a matrix R0 and a vector rdiag based on the struct obsrates,
%   to enable parametrization of the measurement model of the EngEr
%   Covid-19 model. The struct obsrates has the following fields:
%   
%   -states: Cell array where each element define the indices of
%   states which are summed to produce a measurement in the model,
%   i.e. [1,2] produces a measurement of the sum of state#1 and state
%   #2
%
%   -indobs: Cell array where each element define the indices of
%   states which contribute to a measurement n the model through a
%   linear combination.
%
%   -indobspars: Cell array where each element define the weights of
%   the states defined in indobs to determine the measurement, i.e an
%   element [1,2] in indobs and the corresponding element [a,b] in
%   indobspars indicates a measurement of the form a*state#1 +
%   b*state#2.
%
%   -nstate: Number of states in total in the model.
%
%   -R0: Array of additive measurement noise
%
%   -rdiag: Array of multiplicative measurement noise
%
%   -R0_p2: Array of additive measurement noise for after change to
%   experimental sewage method
%
%   -rdiag_p2: Array of multiplicative measurement noise after
%   change to experimental sewage method
%   
%   The output H is a sparse matrix with obsrate.nstate columns where
%   every row represents one measurement. The measurements defined by
%   obsrates.states define the first numel(obsrates.states) rows, the
%   remaing numel(obsrates.indobs) rows are given by obsrates.indobs
%   and obsrates.indobspars. The outpt R0 is a sparse diagonal matrix
%   with diagonal elements defined by obsrates.R0, while rdiag is
%   given directly by obsrates.rdiag.

% J. Evaeus 2022-10-24 - Allow for change in sewage method
% H. Runvik 2020-10-14

% direct measurements (weight = 1)
ii = repelem(1:numel(obsrates.states),cellfun('prodofsize',obsrates.states));
jj = cat(2,obsrates.states{:});
ss = ones(size(jj));

% merge with indirect measurements (weights are given)
if ~isempty(obsrates.indobs)
    ii = [ii repelem(numel(obsrates.states)+(1:numel(obsrates.indobs)), ...
                 cellfun('prodofsize',obsrates.indobs))];
    jj = [jj cat(2,obsrates.indobs{:})];
    ss = [ss cat(2,obsrates.indobspars{:})];
end

% final assemble
H = sparse(ii,jj,ss,max(ii),obsrates.nstate);

% covariance model assuming only diagonal elements
R0 = spdiags(obsrates.R0,0,max(ii),max(ii));
rdiag = obsrates.rdiag;

if isfield(obsrates, 'R0_p2')
  R0_p2 = spdiags(obsrates.R0_p2,0,max(ii),max(ii));
  rdiag_p2 = obsrates.rdiag_p2;
else
  R0_p2 = NaN;
  rdiag_p2 = NaN;
end