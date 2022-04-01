%LAGGEN generates the lag data for the figures and predictions using the
% posterior and Kalman filter.

% snippets of code for copy and paste (for weekly report)
% defaults
posteriordate = '210531'; % change here
ending        = '1_100'; % at what date does the slabs start.
datadate     = [200401 210531] % change end date here

datasource   = 'RU'; % Employs C19 if not Uppsala, then Uppsala Region
useCSSS      = false;

% run all regions
type = 1, lag = 14;
regionList = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
  'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
  'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
  'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
  'Jämtland' 'Västerbotten' 'Norrbotten'};
for reg = [1:21]
  region = regionList{reg}
  FINALRUN = true;
  weekly_prediction
end


% lag7-plot, Uppsala
FINALRUN = true;
type = 2, lag = 7;
region = 'Uppsala'
weekly_prediction


% lag7-plot, Stockholm
FINALRUN = true;
type = 2, lag = 7;
region = 'Stockholm'
weekly_prediction

% lag14-plot, Sweden
FINALRUN = true;
type = 1, lag = 14;
region = 'Sweden'
weekly_prediction

% lag7-plot, Sweden
FINALRUN = true;
type = 2, lag = 7;
region = 'Sweden'
weekly_prediction

disp('*** Done ***');
