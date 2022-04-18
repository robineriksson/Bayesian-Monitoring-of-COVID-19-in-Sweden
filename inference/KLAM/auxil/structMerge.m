function [rates] = structMerge(rates0,rates1,keep0,keep1)
% RATES = structMerge(RATES0,RATES1,KEEP0,KEEP1) merges RATES1 onto RATES0
%   and while keepinging only KEEP0 [rows, cols] from RATES0 and KEEP1
%   [rows, cols] form RATES1.
%
%   The output, RATES, will have size
%   size(RATES0)-KEEP0+size(RATES1)-KEEP1. Where RATES0 and RATES1 are of
%   the same type.
%
%   KEEP0 = {rows, cols}, cell of integer arrays.
%   KEEP1 = {rows, cols}

% R Eriksson 2020-11-11

rates = {};
names = fieldnames(rates0);
names_ = fieldnames(rates1);

% find dynamic names
names_unique = unique(names);
names_unique_ = unique(names_);
hasmeta = find(strcmp(names_unique,'meta'));
hasmeta_ = find(strcmp(names_unique_,'meta'));
if ~isempty(hasmeta_)
  names_unique_(hasmeta_) = [];
end
if ~isempty(hasmeta)
  names_unique(hasmeta) = [];
end

if ~isequal(names_unique,names_unique_)
  error('cannot merge two structs');
end

slabdependent = {};
for i = 1:numel(names_unique)
  name = names_unique{i};
  if size(rates0.(name),1) > 1 ||  size(rates1.(name),1) > 1
    slabdependent{numel(slabdependent)+1} = name;
  end
end

% merge
exceptlist = {'meta'};
for i = 1:numel(names)
  name = names{i};
  if keep0{1} == 0 % merge only static ones
    if any(strcmp(name,slabdependent)) % dynamic
      rates.(name) = rates1.(name)(keep1{1},keep1{2});
    else % static
      rates.(name) = rates0.(name);
    end
  else % merge all
    if ismember(name,fieldnames(rates1)) && ~ismember(name,exceptlist)
      if size(rates0.(name),1) > 1 % dynamic (per slab) | append rates1
        rates.(name) = [rates0.(name)(keep0{1},keep0{2}); rates1.(name)(keep1{1},keep1{2})];
        % *** append rates1 on rates0 | not a below
        %[~, ix] = sort([keep0{1} keep1{1}]);
        %rates.(name) = rates.(name)(ix,:);
      else
        rates.(name) = rates0.(name)(:,keep0{2}); % static (use rates0)
      end
    else
      warning([name ' is missing in rates1'])
    end
  end
end
end
