%LAGGEN generates the lag data for the figures and predictions using the
% posterior and Kalman filter.
%
% All simulation can be generated using
% reg= [1:22], type = 1, laggen, type=2, laggen;
%
% But all data is often not needed, and the specific region can be
% specified by REG and the type of output, 7-day lag ahead (TYPE = 1)
% and filter output with 14 day ahead prediction at data stop (TYPE = 2).
%

% run all regions
if ~exist('reg','var')
    reg = [1:22];
end
if ~exist('type','var')
    type = 1; lag = 14;
end

if type == 1
    lag = 14;
elseif type == 2
    lag = 7;
else
    error('plot type specified is not defined');
end

if ~exist('register','var')
    register='C19';
end

if ~exist('verb','var')
    verb=false;
end

if ~exist('interp','var')
    interp=2;
end

if interp == 1
    ending = '1_100';
elseif interp == 2
    ending ='1_100_update.mat';
else
    error('Only supporting interp 1 or 2')
end

% defaults
posteriordate = '210531'; % change here

datadate     = [200401 210531]; % change end date here
if verb
    datadate
end

useCSSS      = false;
saveall      = true;
noNetwork    = true;


regionList_ = regions(false);
regionList_ = cat(1,regionList_,'Sweden');
for rid = reg
    region = regionList_{rid};
    if verb
        disp(['running: ' region]);
        disp(['with register: ' register]);
        disp(['and datadate: ' num2str(datadate)])
    end
    FINALRUN = true;
    weekly_prediction
end

if verb
    disp('*** Done ***');
end
