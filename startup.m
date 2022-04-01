%STARTUP Startup for C19KALMAN.

% S. Engblom 2021-04-01 (partially new structure)
% S. Engblom 2020-09-16

% link = location of this startup.m
link = mfilename('fullpath');
link = link(1:end-7); % s-t-a-r-t-u-p is 7 chars

% path to folders
addpath(genpath([link 'data/']));
addpath(genpath([link 'inference/']));
addpath(genpath([link 'kalman/']));
addpath(genpath([link 'paper/']));
addpath(genpath([link 'test/']));
addpath(genpath([link 'URDME/']));
addpath(genpath([link 'weekly/']));
% (note: don't add to old/)