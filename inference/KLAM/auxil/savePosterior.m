function [rates, rates100, ratesR0] = savePosterior(thetas, sl, slabs, amparam, ...
                                burnin,jump, perRegion, useCSSS,tosave,...
                                fix,slabtype,smc,id,type)
%RATES = SAVEPOSTERIOR saves the posterior in the correct format.
%   SAVEPOSTERIOR formulates the posterior file that is generated using
%   KLAM and then loaded by other functions, e.g., posteriorenger.m
%
%   *Input*
%   THETAS    - posterior matrix.
%   KL        - (synthetic) likelihood for each slab and sample.
%   SLAB      - The slabs used in the posterior.
%   AMPARAM   - parameters, and others, used in the AM algorithm, as data
%   info.
%   BURNIN    - what burnin length to use when storing the markov chain.
%   JUMP      - For large files, we cannot store all samples, only store
%   every JUMP value.
%   PERREGION - Was the posterior generated using single region data or
%               transports?
%   USECSSS   - if CSSS data was used durring the training, csss tag is
%               then added to the name.
%   TOSAVE    - true/false, if false only return don't save the file.
%
%   FIX       - true/false, if sigma=gammaI or sigma!=gammaI
%
%   SLABTYPE  - {0,1,2,3}, where are the slab cuts placed?
%
%   *Output*
%   RATES    - the struct that is stored in the display message.
%   RATES100 - 100 long version of (thinned) rates.
%   RATESR0  - only R0 from rates.

% R. Eriksson 2021-02-23
if nargin < 8
  burnin=5e4;
  jump=1e2;
  perRegion=false;
  useCSSS=false;
  tosave=true;
  fix=false;
  slabtype=0;
  smc=0;
  id=[];
end

if smc
  rates = mat2struct(thetas(:,:,end)',amparam.ratenames,...
    amparam.hyp, amparam.fix,amparam.nslab);
else
  rates = mat2struct(thetas(:,burnin:jump:end),amparam.ratenames, ...
    amparam.hyp, amparam.fix,amparam.nslab);
  d = round((size(thetas,2)-burnin-1)/100);
  rates100 = mat2struct(thetas(:,burnin:d:end),amparam.ratenames, ...
    amparam.hyp, amparam.fix,amparam.nslab);

  ratesR0 = struct('R0',{rates.R0});
end

%% .meta
% meta information, e.g., hash, Fhash, revision, ...
meta = struct();

% .prorRev
meta.priorRev = amparam.priorRev;

% .revision
meta.rev = datestr(today());

% data .reg/rev/hash
meta.dataReg = amparam.dataReg;
meta.dataRev = amparam.dataRev;
meta.dataHash = amparam.dataHash;

% Other fields that would explain the posterior further.
meta.date = amparam.date;
meta.slabstop = amparam.slabstop;
meta.slabs = slabs;
meta.region = amparam.region;
try
  meta.hypfile = amparam.hypfile;
end


% .hash
cellrates = struct2cell(rates);
cellrates = cat(1,cellrates{:});
cellrates = cellrates(:);
meta.hash = fsetop('check',cellrates);

% .origo
try
  meta.origo = amparam.origo.file;
catch
  ; %do nothing.
end

% meta100, metaR0 - the hash is a bit different for the reduced sizes
meta100 = meta;
cellrates100 = struct2cell(rates100);
cellrates100 = cat(1,cellrates100{:});
cellrates100 = cellrates100(:);
meta100.hash = fsetop('check',cellrates100);


metaR0 = meta;
cellratesR0 = struct2cell(ratesR0);
cellratesR0 = cat(1,cellratesR0{:});
cellratesR0 = cellratesR0(:);
metaR0.hash = fsetop('check',cellratesR0);

% append meta to rates structs
rates.meta = meta;
rates100.meta = meta100;
ratesR0.meta = metaR0;



% likelihood
sl_burn = sl(:,burnin:jump:end);
sl_burn100 = sl(:,burnin:d:end);

filename = mfilename('fullpath');
filename = [filename(1:end-24) 'results/'];
if perRegion
  filename = [filename 'KLAM/perRegion/'];
end

if smc
  filename = [filename 'smc'];
else
  filename = [filename 'slam'];
end
filename = [filename  num2str(amparam.date(end))];

filename = [filename '_' amparam.region '_monthly'];
if useCSSS
  filename = [filename '_csss'];
end
if fix
  filename = [filename '_fix'];
end
switch slabtype
  case 1
    filename = [filename '_1'];
  case 8
    filename = [filename '_8'];
  case 15
    filename = [filename '_15'];
  case 22
    filename = [filename '_22'];
  otherwise
    error('not defined slabtype');
end

if ~isempty(id)
  filename = [filename '_run' num2str(id)];
end

if contains(meta.dataReg,'URDME')
  filename = [filename '_' meta.dataReg];
end

filename100 = [filename '_100.mat'];
filenameR0 = [filename '_R0.mat'];
filename = [filename '.mat'];
if tosave
 if ismember('full', type)
  save(filename, 'rates', 'slabs', 'sl_burn', 'amparam');
  disp(['saved: ' filename]);
 end

 if ismember('100', type)
  rates = rates100;
  save(filename100, 'rates', 'slabs', 'sl_burn100', 'amparam');
  disp(['saved: ' filename100]);
 end

 if ismember('R0', type)
  rates = ratesR0;
  save(filenameR0, 'rates', 'sl_burn', 'slabs', 'amparam');
  disp(['saved: ' filenameR0]);
 end
else % don't save, only return
 disp(['did not save: ' filename]);
  disp(['did not save: ' filename100]);
  disp(['did not save: ' filenameR0]);
end
end