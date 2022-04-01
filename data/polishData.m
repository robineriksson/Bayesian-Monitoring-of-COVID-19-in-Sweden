function data = polishData(data,cumfield,incfield,type)
%DATA = POLISHDATA(DATA,CUMFIELD,INCFIELD,TYPE) polish data.
%   DATA = POLISHDATA(DATA,CUMFIELD,INCFIELD) checks and ensures that
%   there is no negative incidence data. Negative values in the fields
%   INCFIELD are detected and corrected to 0, with values backwards in
%   time corrected in a "reasonable" way, and with CUMFIELD adjusted
%   accordingly.
%
%   DATA follows the format described in LOADDATA. Both INCFIELD and
%   CUMFIELD are cell-vectors of corresponding strings.
%
%   DATA = POLISHDATA(DATA,...,TYPE) additionally replaces unreliable
%   measurements from C19-data with NaNs (TYPE = 0) or with linearly
%   interpolated data from adjacent days (TYPE = 1, the default). Use
%   TYPE = NaN or empty to skip this additional check.
%
%   Example:
%     Data = loadData('C19'); Data0 = Data;
%     Data = polishData(Data,'D','Dinc');
%
%     region = 4; % 4, 6, 10 are the worst (then 3 & 14)
%     d = Data.Dinc(:,region)-Data0.Dinc(:,region);
%     figure(1), clf, stairs(d);
%     figure(2), clf, stairs(cumsum(d,'omitnan'));
%
%   See also LOADDATA.

% The elements that are replaced were determined via manual
% inspection of the C19 data performed 20-11-16 (re-checked 21-02-26).
%
% Data = loadData('C19');
% for k1 = 1:21
%   figure(1), plot(Data.H(:,k1))
%   figure(2), plot(Data.W(:,k1))
%   figure(3), plot(Data.D(:,k1))
%   title(num2str(k1))
%   pause
% end

% S. Engblom 2021-02-26 (Revision, new algorithm)
% S. Engblom 2021-02-25 (merged checkNegInc & removeUnrealiableC19ata)

if nargin < 4
  type = 1;
  if nargin < 3
    incfield = {};
    if nargin < 2
      cumfield = {};
    end
  end
end
if ~iscell(cumfield), cumfield = {cumfield}; end
if ~iscell(incfield), incfield = {incfield}; end
if numel(cumfield) ~= numel(incfield)
  error('Can not match cumulative and incremental fields.');
end

% add a note that data has been polished
if contains(data.rev,'polishedData')
  warning('Data seems to have already been polished once.');
end
data.rev = [data.rev ' | polishData'];

% (1) find negative values in incfield, set them to 0 and adjust
% smoothly Ndays backwards in time; repeat recursively if need be

% parameter Ndays: time window backwards within which the fix will
% (recusively) take place
Ndays = 14;
for ix = 1:numel(incfield)
  % for convenience
  Iinc = data.(incfield{ix});

  % "credit and debit": separate into (non-NaN) positive and negativa parts
  IincNAN = isnan(Iinc); % (will be put back later)
  Iinc(IincNAN) = 0;
  Iincpos = max(Iinc,0);
  Iincneg = -min(Iinc,0);
  % loop over columns containing negative increments
  cols = find(any(Iincneg,1));
  for jj = cols
    % loop *backwards* in each column
    ii = find(any(Iincneg(:,jj),2),1,'last');
    while ~isempty(ii)
      % Ndays "smoothed out" corrector
      corr = repmat(floor(Iincneg(ii,jj)/Ndays),Ndays,1); % equal part...
      rem = mod(Iincneg(ii,jj),Ndays); % ...remainder to be distributed...
      corr = corr+repelem([1 0],[rem Ndays-rem])'; % ...in unequal parts
      % find the last correction
      iend = find(corr,1,'last');
      corr = corr(1:iend);
      if ii-iend <= 0
        % attempt to correct before the beginning of the array
        error('Correction could not be performed within window.');
      end
      % remove the offender and find the first meaningful index to correct
      Iincneg(ii,jj) = 0;
      istart = find(Iincpos(ii-1:-1:1,jj),1,'first');
      % correct in window backwards
      Iincpos(ii-istart:-1:ii-istart-iend+1,jj) = ...
          Iincpos(ii-istart:-1:ii-istart-iend+1,jj)-corr;
      % adjust credit and debit accordingly
      Iincneg = Iincneg-min(Iincpos,0); % more negatives might creep in
      Iincpos = max(Iincpos,0); % should remain positive at all times
      ii = find(any(Iincneg(:,jj),2),1,'last');
    end
  end
  % assemble final
  Iinc = Iincpos-Iincneg;
  if any(Iinc(IincNAN),'all')
    error('Correction within window could not be performed before NaNs.');
  end
  Iinc(IincNAN) = NaN;

  % restore result to fields
  data.(incfield{ix}) = Iinc;
  data.(cumfield{ix}) = cumsum(data.(incfield{ix}),1,'omitnan');
end

% (2) fix of suspicious elements (only for C19/RU)
if isempty(type) || isnan(type) || ...
      ~any(strcmp(data.reg,{'C19' 'RU'}))
  return;
end

% manually found suspicious H-measurements
NaNidx_H = data.H < 0; % Find negative values.
NaNidx_H(60:62,1) = true;
NaNidx_H(47:48,2) = true;
NaNidx_H(47:53,14) = true;
NaNidx_H(55,18) = true;
NaNidx_H(78,19) = true;
NaNidx_H(47:48,20) = true;
NaNidx_H(55,20) = true;
NaNidx_H(48,21) = true;

NaNidx_W = NaNidx_H; % H and W incorrect at same days...
NaNidx_W(47:48,14) = false; % ...except for region 14 which is odd
%NaNidx_W(60:62,1) = false;

switch type
 case 0
  data.H(NaNidx_H) = NaN;
  data.W(NaNidx_W) = NaN;
 case 1
  xgrid = (1:numel(data.date))';
  for j = 1:numel(data.regions)
    % H: create grid from x-values, removing all NaN-coordinates
    xh = xgrid;
    xh(NaNidx_H(:,j)) = [];
    % interpolate the NaN-coordinates
    data.H(NaNidx_H(:,j),j) = interp1(xh,data.H(~NaNidx_H(:,j),j), ...
                                      xgrid(NaNidx_H(:,j)));
    % ditto W
    xw = xgrid;
    xw(NaNidx_W(:,j)) = [];
    data.W(NaNidx_W(:,j),j) = interp1(xw,data.W(~NaNidx_W(:,j),j), ...
                                      xgrid(NaNidx_W(:,j)));
  end
end
