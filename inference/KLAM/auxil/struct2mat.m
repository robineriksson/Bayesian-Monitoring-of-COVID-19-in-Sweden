function [x0, ratenames, meta] = struct2mat(rates)
%[X0] = STRUCT2MAT(RATES) convert a struct RATES into a matrix X0
%   RATES is a struct, and to be able to use the values efficiently in the
%   Adaptive Metropolis algorithm, we need to be able to convert the full
%   rates struct into a matrix.
%
%   *input*
%   RATES - rates struct, see priorenger/posteriorenger
%
%   *output*
%   X0 - appended rates in matrix form.
%   RATENAMES - fieldnames of rates.
%   META - struct that hold the meta variables of the rates structure:
%   {hash Fhash rev}.

% R. Eriksson 210522

fnames = fieldnames(rates);
x0 = [];

% we can exlude the F* cases as these can be loaded in and formated again
% later.
exclude = {'F0ave' 'F1ave' 'F2ave' 'F2dave' 'F3ave' 'F3dave' 'F4ave' ...
           'F2scale' 'F2dscale' 'F3scale' 'F3dscale' ...
           'F1' ...
           'SIR_MORT' ...
           'rev' 'reg' 'chunks' 'hash'...
           'dataReg' 'dataRev' 'dataHash' ...
           'Fhash' 'priorRev' 'priorHash' ...
           'A2I' ...
           'meta'};
% which ones has size which changes dependent on the the number of slabs?
slabdependent = {};
for i = 1:numel(fnames)
  name = fnames{i};
  if size(rates.(name),1) > 1
    slabdependent{numel(slabdependent)+1} = name;
  end
end

include = fnames(~fsetop('ismember',fnames, exclude));

% append the parameters.
for i = 1:numel(include)
  name = include{i};
  X = rates.(name);
  x0 = [x0; X];
end

ratenames = cell(size(x0,1),1);

% keep track on row belongs to what rates.
times = size(rates.lambda,1)-1;
shortlist = find(fsetop('ismember',include, slabdependent));


name_id = 1;
slot_id = 1;
while slot_id <= size(x0,1)
  name = include(name_id);
  if any(name_id == shortlist)
    ratenames(slot_id:(slot_id+times)) = name;
    slot_id = slot_id+times+1;
  else
    ratenames(slot_id) = name;
    slot_id = slot_id+1;
  end
  name_id = name_id + 1;
end

if any(strcmp(fieldnames(rates), 'hash'))
  meta = struct;
  try
    meta.hash = rates.hash;
  catch
    meta.hash = [];
  end

  try
    meta.Fhash = rates.Fhash;
  catch
    meta.Fhash = [];
  end

  try
    meta.rev = rates.rev;
  catch
    meta.rev = [];
  end
else
 meta = rates.meta;
end

end
