%LAGGEN generates the lag data for the figures and predictions using the
% posterior and Kalman filter.

% run all regions
if ~exist('reg','var')
    reg = [1:22]
end
if ~exist('type','var')
    type = 1, lag = 14;
end

if type == 1
    lag = 14;
elseif type == 2
    lag = 7;
else
    error('plot type specified is not defined');
end

if ~exist('register','var')
    register='C19'
end
% defaults
posteriordate = '210531'; % change here
ending        = '1_100'; % at what date does the slabs start.
datadate     = [200401 210531] % change end date here


useCSSS      = false;
saveall      = true;
noNetwork    = true;


regionList_ = regions(false);
regionList_ = cat(1,regionList_,'Sweden');
for rid = reg
  region = regionList_{rid}
  FINALRUN = true;
  weekly_prediction
end

disp('*** Done ***');
