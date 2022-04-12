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

% defaults
posteriordate = '210531'; % change here
ending        = '1_100'; % at what date does the slabs start.
datadate     = [200401 210531] % change end date here

datasource   = 'C19';
useCSSS      = false;

saveall      = true;

noNetwork    = true;


regionList = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
  'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
  'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
  'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
  'Jämtland' 'Västerbotten' 'Norrbotten' 'Sweden'};
for rid = reg
  region = regionList{rid}
  FINALRUN = true;
  weekly_prediction
end

disp('*** Done ***');
