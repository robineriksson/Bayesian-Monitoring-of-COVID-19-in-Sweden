function state0 = getInitState(useinit, date, region)
%GETINIT returns the initialization state generated from create_C19initfilt
%   STATE0 = GETINITSTATE(true,DATE,REGION) finds the file for a specific
%   DATE, and REGION, and returns a struct in the desired format.
%
%   If USEINIT = false, then the function returns an empty variable. This
%   functionality is included to reduce the size of the runnable scripts.

% R. Eriksson 210522

if useinit
  % folder location
  folder = mfilename('fullpath');
  folder = [folder(1:end-23) 'results/'];
  
  % construct struct
  state0 = struct();
  state0.date = date;
  % define file
  state0.file = ['slam' num2str(state0.date) '_' region '_monthly_1'];
  loadfile = [folder 'state/perRegion/' state0.file '.mat'];
  X0 = load(loadfile,'Y0','Y0Cov');%,'Ydata','date','useCSSS','datasource');
  
  state0.Y0 = X0.Y0;
  state0.Y0Cov = X0.Y0Cov;
else
  state0 = [];
end