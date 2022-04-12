function p = postpath
%POSTPATH Returns the base directory of the posterior distributions.

% S. Engblom 2021-05-27

p = mfilename('fullpath');
p(end-numel(mfilename)+1:end) = [];