%GEN_CREATE_INITC19FILT cals CREATE_INITC19FILT to generates a savepoint
%  (state posterior including hidden states) that allows one to start from
%  a later state, i.e., instead of restarting from data0 this creates a new
%  origo (starting point).
 
% R. Eriksson 05-12-21


% identify the posterior file
posteriordate = '210331';
shift        = '1'; % at what date does the slabs start.
useCSSS      = false; % use CSSS for I(+A) compartment?
use100       = true;
% simulation stop date
date         = 210331; % end date here



% folder for posteriors
abspath = mfilename('fullpath');
prefix = [abspath(1:end-35) 'inference/results/SLAM/'];

regionList = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
 'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
 'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
 'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
 'Jämtland' 'Västerbotten' 'Norrbotten'};

% main loop, iterate through all the regions, and save the files under
% /inference/results/state/
for k = 1:numel(regionList)
 region = regionList{k};
 
 % construct the posterior file name
 posterior = ['perRegion/slam' posteriordate '_' region '_monthly'];
 if useCSSS
  posterior = [posterior '_csss'];
 end
 posterior = [posterior '_' shift];
 if use100
   posterior = [posterior '_100'];
 end
 posterior = [posterior '.mat'];
 
 % generate data
 [Y0, Y0Cov, meta] = create_initC19filt([],[prefix posterior],date);
 
 % save to file
 savefile = [abspath(1:end-35) 'inference/results/state/' ...
             posterior(1:end-8) '.mat'];
 save(savefile,'Y0','Y0Cov','meta');
 disp(['saved: ' savefile]);
end
