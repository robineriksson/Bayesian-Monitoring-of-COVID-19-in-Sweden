function umod = covid19enger2
%COVID19ENGER2 Scenario model file for the EngEr Covid-19 model.
%   UMOD = COVID19ENGER2 returns an URDME model UMOD for the EngEr
%   Covid-19 model in 21 counties.
%
%   The modification here divides the D compartment into D_I, D_H, and
%   D_W.  Additionally, H -> W <-> H2.
%
%   See also COVID19ENGER.

% S. Engblom 2021-04-17 (Major revision)
% S. Engblom 2020-06-26 (Revision)
% S. Engblom 2020-04-09 (Revision)
% S. Engblom 2020-03-22

% the variable phi can be either understood as a state variable 'phi'
% or as a (local) rate 'PHI'

% rates: most are gdata, but beta and PHI are local
rates.beta = 'ldata';
rates.PHI = 'ldata';
% fractions:
rates.F0 = 'gdata';  % E --> I
rates.F1 = 'gdata';  % A --> I
rates.F2 = 'gdata';  % I --> H
rates.F2d = 'gdata'; % I --> D
rates.F3 = 'gdata';  % H --> W
rates.F3d = 'gdata'; % H --> D
rates.F4 = 'gdata';  % W --> D
% actual rates:
rates.sigma = 'gdata';  % E -->
rates.gammaI = 'gdata'; % I -->
rates.gammaA = 'gdata'; % A -->
rates.gammaH = 'gdata'; % H -->
rates.gammaW = 'gdata'; % W -->

% species, transitions and rates
species = {'I' 'A' 'E' 'phi' 'H' 'W' 'D_I' 'R' 'H2' 'D_H' 'D_W'};
r = cell(1,16);
r{1} = '@ > beta*PHI > E'; % PHI (double) or phi (integer)
r{2} = 'E > sigma*F0*E > I';
r{3} = 'E > sigma*(1-F0)*E > A';
r{4} = 'A > gammaA*(1-F1)*A > R';
r{5} = 'A > gammaA*F1*A > I';
r{6} = 'I > gammaI*F2*I > H';
r{7} = 'I > gammaI*F2d*I > D_I';
r{8} = 'I > gammaI*(1-F2-F2d)*I > R';
r{9} = 'H > gammaH*F3*H > W';
r{10} = 'H > gammaH*F3d*H > D_H';
r{11} = 'H > gammaH*(1-F3-F3d)*H > R';
r{12} = 'W > gammaW*F4*W > D_W';
r{13} = 'W > gammaW*(1-F4)*W > H2';
r{14} = 'H2 > gammaH*F3*H2 > W';
r{15} = 'H2 > gammaH*F3d*H2 > D_W';
r{16} = 'H2 > gammaH*(1-F3-F3d)*H2 > R';

% sort it out
umod = rparse([],r,species,rates,'covid19enger.c');

% add demographic information
load Ncounties
umod.private.POP = sum(N,1);
