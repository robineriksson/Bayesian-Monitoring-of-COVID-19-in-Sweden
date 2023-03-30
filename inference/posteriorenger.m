function posterior = posteriorenger(Nreplicas,file,weights,perFile)
%POSTERIORENGER Posterior parameter distribution for the EngEr Covid-19 model.
%   POSTERIOR = POSTERIORENGER(NREPLICAS,FILE) returns NREPLICAS
%   random samples from the posterior FILE. Sampling of the NREPLICAS
%   parameters is a straightforward sampling with replacement, see
%   RANDSAMPLE.
%
%   FILE is the name of the posterior .mat file. The naming should
%   preferably follow the pattern 'abcYYMMDD_XX.mat', where XX
%   corresponds to [low/mid/high, low/mid/high] and is represented by
%   [1/2/3 1/2/3]. The file should preferably be stored in
%   /inference/results/.
%
%   POSTERIOR = POSTERIORENGER(inf,FILE) returns the mean of the
%   posterior.
%
%   POSTERIOR = POSTERIORENGER([],FILE) returns all available samples.
%
%   POSTERIOR = POSTERIORENGER(...,{F1 F2 ...},[W1 W2 ...]) produces a
%   single sample from all files F1, F2, ... weighted according to the
%   weights W1, W2, ...
%
%   The fourth input PERFILE {false} triggers, if true, per file R0
%   (and beta), i.e., size(POSTERIOR.R0,3) == numel(FILE)
%
%   In all the above syntaxes, FILE may also be a rate struct or cell
%   vector of such. This is mainly intended for testing purposes but
%   may occasionally be useful in other situations.
%
%   Example:
%     % basic use: 100 samples
%     file = 'slam210418_22_Dalarna_monthly_15';
%     P = posteriorenger(100,file);
%
%     % 'combined' syntax with weights based on county population
%     files = {'slam210418_22_Dalarna_monthly_15' ...
%              'slam210418_22_Halland_monthly_15'};
%     load Ncounties;
%     W = sum(N,1);
%     P = posteriorenger(100,files,W([16 11])); % 16: Dalarna, 11: Halland
%     % weighted mean:
%     M = posteriorenger(inf,files,W([16 11])); % 16: Dalarna, 11: Halland
%
%   See also PRIORENGER, RANDI, RANDSAMPLE.

% S. Engblom 2021-04-23 (Revision, new output)
% R. Eriksson 2020-10-12 (New assumption on saving format)
% R. Eriksson 2020-09-30

if nargin < 4
  perFile = false;
end
% 'combined' syntax
if nargin > 2 % == 3
  % 'all'
  if isempty(Nreplicas)
    % simply concatenate all available samples
    posterior = posteriorenger([],file{1});
    fnames = fieldnames(posterior);
    metanames = fieldnames(posterior.meta);
    % ...then continue with the rest
    for i = 2:numel(file)
      p = posteriorenger([],file{i});
      posterior = l_catstruct(posterior,p,fnames,metanames);
    end
    % (note: this object cannot be compared to a full urn as the samples
    % are not weighted correctly)
    return;
  end

  % weights per file
  npost = numel(weights);
  weights = [0 cumsum(weights)];
  weights = weights/weights(end);

  if Nreplicas == inf
    % weighted means: means first...
    posterior = posteriorenger(inf,file{1});
    fnames = fieldnames(posterior);
    metanames = fieldnames(posterior.meta);
    for i = 2:numel(file)
      p = posteriorenger(inf,file{i});
      posterior = l_catstruct(posterior,p,fnames,metanames,perFile);
    end
    % ...then weight them together
    if perFile % Don't compute the average of the R0 values.
      exclude = {'meta' 'R0'};
    else
      exclude = {'meta'};
    end

    for j = 1:numel(fnames)
      if ~any(strcmp(fnames{j},exclude))
        posterior.(fnames{j}) = posterior.(fnames{j})*diff(weights(:));
      end
    end
    if perFile % region based R0 with weighted mean of the other rates.
      posterior.beta = getC19beta(posterior,posterior.R0,posterior.meta.interp);
    end
    return;
  end

  % otherwise: main syntax, sample from several urns according to weight

  % sample the number of samples per file
  count = histc(rand(1,Nreplicas),weights);
  count = count(1:end-1);

  % sample from the first...
  posterior = posteriorenger(count(1),file{1});
  fnames = fieldnames(posterior);
  metanames = fieldnames(posterior.meta);
  % ...then continue with the rest
  for i = 2:numel(file)
    p = posteriorenger(count(i),file{i});
    posterior = l_catstruct(posterior,p,fnames,metanames);
  end
  return;
end

% load all samples
if ischar(file)
  load(file,'rates');
elseif isstruct(file)
  rates = file;
else
  error('Unknown posterior file format.');
end

% all/random sampling/compute average
fnames = fieldnames(rates);

% meta field is not sampled:
[~,fsample] = fsetop('setdiff',fnames,{'meta'});
posterior = struct();

if isempty(Nreplicas)
  % 'all'-syntax
    for i = fsample
        posterior.(fnames{i}) = rates.(fnames{i});
    end
elseif Nreplicas < inf
  % sampling Nreplicas *with* replacement
%  idout = randsample(1:size(rates.(fnames{fsample(1)}),2),Nreplicas,'true');
  % faster:
  idout = randi(size(rates.(fnames{fsample(1)}),2),1,Nreplicas);

  for i = fsample
      posterior.(fnames{i}) = rates.(fnames{i})(:,idout);
  end
else
  % the average of the posterior
  for i = fsample
    posterior.(fnames{i}) = mean(rates.(fnames{i}),2);
  end
end

% copy .meta
posterior.meta = rates.meta;

% if there is any origo file -> merge!
if isfield('origo',rates.meta)
  origofile = rates.meta.origo;
  origo = posteriorenger([],[origofile '_100']);

  % merge rates
  fnames = fieldnames(origo);

  % when mixing, we need to know how many of each rate file to use.
  w0 = numel(origo.meta.date); % weight origo
  w = numel(posterior.meta.date); % weight posterior
  ws = w0+w; % normalize
  w0 = w0/ws; w = w/ws;
  N = size(posterior.(fnames{fsample(1)}),2);
  for k = 1:numel(fnames)
    if ~any(strcmp(fnames{k},metafields))
      if size(origo.(fnames{k}),1) > 1 % dynamic | just append new on origo
        idout = randi(size(origo.(fnames{k}),2),1,N);
        posterior.(fnames{k}) = [origo.(fnames{k})(:,idout); ...
                                 posterior.(fnames{k})];
      else % fixed | merge based on datasize.
        % merge based in data size.
        idout0 = randi(size(origo.(fnames{k}),2),1,round(N*w0));
        idout = randi(size(posterior.(fnames{k}),2),1,round(N*w));
        if numel(idout0) + numel(idout) ~= N
          error('rounding problem. output not equal to N');
        end
        posterior.(fnames{k}) = [origo.(fnames{k})(idout0), ...
                                 posterior.(fnames{k})(idout)];
      end
    end
  end

  % merge meta
  % .date
  posterior.meta.date = [origo.meta.date; posterior.meta.date];
  % .slabs
  posterior.meta.slabs = [origo.meta.slabs, ...
           origo.meta.slabs(end)+posterior.meta.slabs];
  % .slabstop
  posterior.meta.slabstop = [origo.meta.slabstop(1:end-1); ...
               origo.meta.slabstop(end)+1;
               origo.meta.slabstop(end)+posterior.meta.slabstop];
  % .hash
  posterior.meta.hash = fsetop('check',...
                               [origo.meta.hash; posterior.meta.hash]);

end

%-----------------------------------------------------------------------
function p = l_catstruct(p1,p2,fnames,metanames,page)
%L_CATSTRUCT Concatenates structs.
%   P = L_CATSTRUCT(P1,P2,FNAMES,METANAMES) concatenates all fields
%   P1.(FNAMES{i}) and P2.(FNAMES{i}) and similarly with the fields
%   P1.meta.(METANAMES{i}) and P2.meta.(METANAMES{i}). If PAGE is given
%   {false}, the R0 field is put on another page size(P.R0) =
%   [nslabs,nsamples,pages]

if nargin < 5
  page = false;
end

% concatenate to the right
for j = 1:numel(fnames)
  if strcmp(fnames{j},'R0')
    if page % add to new page
      p1.(fnames{j}) = cat(3,p1.(fnames{j}), p2.(fnames{j}));
    else % add regularly
      p1.(fnames{j}) = cat(2,p1.(fnames{j}), p2.(fnames{j}));
    end
  elseif ~strcmp(fnames{j},'meta')
      p1.(fnames{j}) = [p1.(fnames{j}) p2.(fnames{j})];
  end
end

% merge .meta
for j = 1:numel(metanames)
  if ischar(p1.meta.(metanames{j}))
    % strings: this looks a bit better:
    p1.meta.(metanames{j}) = ...
        [p1.meta.(metanames{j}) ' + ' p2.meta.(metanames{j})];
  else
    if ~fsetop('ismember',metanames(j),{'date' 'slabstop' 'slabs'})
      % fields can be concatenated:
      p1.meta.(metanames{j}) = ...
          [p1.meta.(metanames{j}) p2.meta.(metanames{j})];
    else
      % fields not concatenated but checked:
      if ~isequal(p1.meta.(metanames{j}),p2.meta.(metanames{j}))
        error(sprintf('Field .%d mismatch. Cannot concatenate.',metanames{j}));
      end
    end
  end
end
p = p1;

%-----------------------------------------------------------------------
