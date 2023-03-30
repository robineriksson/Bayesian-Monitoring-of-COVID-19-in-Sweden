function Data = loadData(regName,quiet)
% LOADDATA(REGNAME) loads C19-data into struct.
%   DATA = LOADDATA(REGNAME) returns a struct DATA with data from the
%   specified source.
%
%   On return, DATA is a struct with fields [model_compartment,
%   regions, date], where,
%     * model_compartment can be {Icum, H, W, Wcum, D, Dinc}
%     * regions are either {'Stockholm', 'Uppsala', ...} or {'Sweden'}
%     * date is of a vector of numeric type: YYMMDD.
%
%   Name             Description
%   ----------------------------------------------------------------
%   'C19'            C19 scraped registry
%   'RU'             C19 + closed source 'Region Uppsala'
%   'FHM'            FHM by region, public registry
%   'FHM_Iinc'       FHM daily incidence confirmed cases
%   'FHM_swe'        FHM daily, Sweden, public registry
%   'FHM2'           FHM, closed source
%   'SS'             Socialstyrelsen, patients in hospital
%   'SS_swe_mort'    Socialstyrelsen, daily dead, Sweden
%   'SIR'            Swedish intensive care registry
%   'SIR_swe_weekly' in ICU (aggregated, Sweden, weekly)
%   'SIR_swe_mort'   30 day Mortality rate @ ICU, Sweden
%   'SIR_swe'        Aggregated reported daily ICU in Sweden
%   'SimInf'         SimInf simulator generated data
%   'CSSS'           Covid symptom study, county estimate
%   'CSSS_public'    Covid symptom study (public github), county estimate
%   'CSSS_swe'       Covid symptom study, Sweden estimate (Uppsala adjusted)
%   'FHM_mutate'     Estimated % of suspected B.1.1.7 variants
%   'FHM_R'          Estimated R per region
%   'WW'             Reported levels of mRNA in wastewater reported per county.
%   'FHM_vac'        FHM, # of vaccinated, 1st and 2nd dose, per region
%   'OWID_swe'       Number of tests administered in Sweden total
%   'FHM_test'       Number of tests per 1000 per region
%   'FHM_prev_swe'   Estimated prevalence in Sweden given by FHM.
%
%    For further information, see the file source.csv.

%    Hidden support for ' -update' appended to all names above. This
%    will update the hash automatically if a mismatch is
%    found. Intended for internal tests and convenience mainly.
%
%    Hidden support also for a second argument QUIET. Set to true to
%    avoid some information being displayed. For testing purposes.
%
%    The logic with hashes is complicated but should appear as rather
%    natural and careful. It can be summarized as follows:
%
%    (§1) first try to load from .mat-file directly;
%         return successfully if doable.
%    (§2) otherwise we need to do a full and much slower load;
%    (§2a) no old hash was found so a new hash is generated;
%         update files and return successfully.
%    (§2b) an old hash was found;
%    (case 1) hash mismatch & ~'-update': will NOT update for fast load
%    (case 2) hash mismatch & '-update': WILL update
%    (case 3) hash match: WILL update
%         THEN return successfully.

% R. Eriksson 2021-02-22 (Fast loading through hash and .mat)
% S. Engblom 2020-10-28 (Revision, more checks)
% R. Eriksson 2020-09-17

if nargin < 2
  quiet = false;
  if nargin < 1
    regName = 'C19';
  end
end

% update syntax?
upd = ' -update';
update = strncmp(regName(end:-1:1),upd(end:-1:1),numel(upd));
if update
  regName = regName(1:end-numel(upd));
end

% parse .csv-file
regName_ = regName; % (used to store hash)
switch regName
  case 'C19'
    fileName = 'C19/C19.csv';
  case 'RU'
    fileName = 'RU/C19RU.csv';
  case 'FHM'
    fileName = 'FHM/FHM_weekly.csv';
  case 'FHM_Iinc'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_daily.csv';
  case 'FHM_swe'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_daily_swe.csv';
  case 'FHM2'
    regName_ = 'FHM';
    fileName = 'FHM/FHM.csv';
  case 'SS'
    fileName = 'SS/SS.csv';
  case 'SS_swe_mort'
    regName_= 'SS';
    fileName = 'SS/SS_swe_mort.csv';
  case 'SIR'
    fileName = 'SIR/SIR.csv';
  case 'SIR_swe_weekly'
    regName_ = 'SIR';
    fileName = 'SIR/SIR_weekly.csv';
  case 'SIR_swe_mort'
    regName_ = 'SIR';
    fileName = 'SIR/SIR_swe_mort.csv';
  case 'SIR_swe'
    regName_ = 'SIR';
    fileName = 'SIR/SIR_swe.csv';
  case 'CSSS'
    regName_ = 'CSSS';
    fileName = 'CSSS/CSSS.csv';
  case 'CSSS_public'
    regName_ = 'CSSS';
    fileName = 'CSSS/CSSS_public.csv';
  case 'CSSS_swe'
    regName_ = 'CSSS';
    fileName = 'CSSS/CSSS_swe.csv';
  case 'SimInf'
    regName = 'SimInf';
    fileName = 'SimInf/SimInf.csv';
  case 'FHM_mutate'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_mutate.csv';
  case 'FHM_R'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_R.csv';
  case 'WW'
    regName_ = 'UWW';
    fileName = 'UWW/UWW.csv';
  case 'FHM_vac'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_vac.csv';
  case 'OWID_swe'
    regName_ = 'OWID';
    fileName = 'OWID/OWID.csv';
  case 'FHM_test'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_test.csv';
  case 'FHM_prev_swe'
    regName_ = 'FHM';
    fileName = 'FHM/FHM_prev_swe.csv';
  case 'test1'
    fileName = 'test1.csv';
  case 'test2'
    regName_ = '';
    fileName = 'test2.csv';

  otherwise
      error('Unknown source.');
end

% read metainformation from source.csv
fullpath = mfilename('fullpath');
sourcepath = [fullpath(1:end-8) 'sources/sources.csv'];
revTable = importdata(sourcepath);
[~,col] = fsetop('ismember',{'file' 'source' 'date' 'hash'}, ...
                 revTable.textdata(1,:));
if any(col == 0)
  error('Could not parse source.csv.');
end
filecol = col(1);
datecol = col(3);
hashcol = col(4);
row = find(fsetop('ismember',revTable.textdata(:,filecol),{fileName}));
if numel(row) == 0
    error(sprintf('Source = %s, is not included in source table.', regName));
end
% metainformation
revData = revTable.textdata(row,datecol);
revData = datestr(revData{:});
% hash is the only data column
revHash = uint32(revTable.data(row-1)); % row-1 because no header.

try
  % (§1) try to load end-result directly
  datapath = [fullpath(1:end-8) 'sources/' regName_ '/' regName 'data.mat'];
  load(datapath,'Data');

  if strcmp(Data.rev,revData) && ...
        strcmp(Data.reg,regName) && ...
        Data.hash == revHash
    % fast load successful:
    return;
    % NOTE: in this case we simply return the Data found in data.mat. The
    % "full hash" in hash.mat is thus not considered.
  end
  throw;
catch
  % (§2) considerably slower...
  l_disp(quiet,['New load of data for register ''' regName '''...']);
  csvpath = [fullpath(1:end-8) 'sources/' fileName];
  try
      p = readcell(csvpath);
      l_disp(quiet,'   ...done. Will load faster next time.');
  catch
      error(sprintf('Cannot read source = %s.', regName));
  end
end

% prepare output struct
Data = struct();

% read headings
Data.regions = p(1,2:end-1);

if length(Data.regions) > 1
% $$$   % generic ordering and naming of the regions
% $$$   load Ncounties Nlan
% $$$   % remove ' län' and 's län' as appropriate
% $$$   for j = 1:numel(Nlan)
% $$$     Nlan{j}(end-3:end) = '';
% $$$     if Nlan{j}(end) == 's'
% $$$       Nlan{j}(end) = '';
% $$$     end
% $$$   end
  % for speed, here is the end-result of the above construct:
  Nlan = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
          'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
          'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
          'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
          'Jämtland' 'Västerbotten' 'Norrbotten'}';

  % replace underscores with spaces
  for j = 1:numel(Data.regions)
    Data.regions{j}(Data.regions{j} == '_') = ' ';
  end

  % work thru the mismatches
  if numel(Nlan) ~= numel(Data.regions)
    error('Missing counties?');
  end
  y = ismember(Data.regions,Nlan);
  for j = find(~y)
    % these mismatches are considered known:
    switch Data.regions{j}
     case 'Sörmland', Data.regions{j} = 'Södermanland';
     case 'Jämtland Härjedalen', Data.regions{j} = 'Jämtland';
     case 'Västragötaland', Data.regions{j} = 'Västra Götaland';
     case 'Jämtland.Härjedalen', Data.regions{j} = 'Jämtland';
     case 'Västra.Götaland', Data.regions{j} = 'Västra Götaland';
     otherwise,
      error([ Data.regions{j} ' is an unknown county.']);
    end
  end
  % deduce the permutation
  [y,yorder] = ismember(Nlan',Data.regions);
  if ~all(y) || any(sort(yorder)-(1:numel(Nlan)))
    error('Cannot parse naming convention of counties.');
  end
  % reorder accordingly
  Data.regions = Data.regions(yorder);
  if ~all(strcmp(Nlan',Data.regions))
    error('Failed in sorting counties.');
  end
else
  yorder = 1;
end

% 1st column is date, final column is compartment
mat_data = p(2:end,2:end-1);
mat_data = mat_data(:,yorder); % correct the ordering of regions.
% Matlab doesn't recoqnize NA as NaN, replace NA with NaN
tf = cellfun(@(x)strcmp(x,'NA'),mat_data);
mat_data(tf) = {NaN};
compartments = cell2mat(mat_data);
% divide up the different model compartments
model_compartments = p(2:end,end);
u_model_compartments = unique(model_compartments);
for i = 1:length(u_model_compartments)
  comp = u_model_compartments(i);
  Data.(comp{1}) = compartments(strcmp(comp, model_compartments), :);
end

% switch to numerical yymmdd-format
tspan = p(2:end,1);
yy = cellfun(@year,tspan);
mm = cellfun(@month,tspan);
dd = cellfun(@day,tspan);
date = (yy-2000)*10000+mm*100+dd;

% repetitive form of date vector, check that we understand it...
date = reshape(date,[],length(u_model_compartments));
if ~all(date(:,1) == date,'all')
  error('Cannot parse date vector.');
end
% ...then pick the first one
Data.date = date(:,1);

% compute and compare hash of output struct
hash = l_hash(Data);

% metainformation
Data.reg = regName;
Data.rev = revData;

% create hash for all of data (recorded in source.csv)
c = struct2cell(hash);
Data.hash = fsetop('check',cat(2,c{:})');

% for some reason, load needs more detailed path than readcell:
hashpath = [fullpath(1:end-8) 'sources/' regName_ '/' regName 'hash.mat'];
if ~exist(hashpath,'file')
  % (§2a) full load of data, no old hash found so new hash generated
  warning('loadData:hash',sprintf('Could not load hash for register ''%s''.',regName));
  disp('Saving new hash, saving data, and updating sources.csv...');
  hash0 = hash;
  save(hashpath,'hash0');

  % save to file for quick load next time
  save(datapath,'Data');

  % register new hash in source.csv to trigger load next time
  revTable.data(row-1) = Data.hash;
  revTable.textdata(2:end,hashcol) = num2cell(uint32(revTable.data));
  writecell(revTable.textdata,sourcepath);

  % NOTE: in this case we update the files data.mat and source.csv for
  % fast load!
  disp('   ...done. Used from now on.');
else
  % (§2b) full load of data, old hash found:
  load(hashpath,'hash0');
  if ~l_comparehash(hash0,hash);
    if ~update
      % (case 1) hash mismatch & ~update: will NOT update
      disp(sprintf('Hash mismatch for register ''%s''.',regName));
      disp('To generate a new hash, simply remove the old one, e.g.,');
      disp(['  !rm -i ' hashpath]);
      disp('and run anew.');
      disp('Will NOT be able to load fast before this is done!');

      % NOTE: so in this case we will NOT update neither the data.mat nor
      % the source.csv! For this to happen you should remove the old
      % hash.mat.
    else
      % (case 2) hash mismatch & update: WILL update
      l_disp(quiet,sprintf('Hash mismatch for register ''%s''.',regName));
      l_disp(quiet,'Saving new hash, saving data, and updating sources.csv...');
      hash0 = hash;
      save(hashpath,'hash0');
      save(datapath,'Data');
      revTable.data(row-1) = Data.hash;
      revTable.textdata(2:end,hashcol) = num2cell(uint32(revTable.data));
      writecell(revTable.textdata,sourcepath);
      l_disp(quiet,'   ...done. Used from now on.');
    end
  else
    % (case 3) hash match: WILL update
    l_disp(quiet,sprintf('Hash MATCH for register ''%s''.',regName));
    l_disp(quiet,'This indicates missing .mat-file for fast load.');
    l_disp(quiet,'Saving data, and updating sources.csv...');

    % save to file for quick load next time
    save(datapath,'Data');

    % register new hash in source.csv to trigger load next time
    revTable.data(row-1) = Data.hash;
    revTable.textdata(2:end,hashcol) = num2cell(uint32(revTable.data));
    writecell(revTable.textdata,sourcepath);

    % NOTE: in this special case (hash from .mat-file match) we update the
    % files data.mat and source.csv for fast load!
    l_disp(quiet,'   ...done. Used from now on.');
  end
end

%-----------------------------------------------------------------------
function hash = l_hash(Data)
%L_HASH Computes a hash structure from data structure.

% region naming/ordering
hash.regions = fsetop('check',fsetop('check',Data.regions)');

% perform hash month-wise
[months,~,im] = fsetop('unique',floor(Data.date'/100));
hash.nmonths = numel(months);

fnames = fieldnames(Data);
for i = 1:numel(fnames)
  if strcmp(fnames{i},'regions')
    continue;
  end
  % real data field
  hash.(fnames{i}) = zeros(1,numel(months),'uint32');
  for j = 1:numel(months)
    hash.(fnames{i})(j) = ...
        fsetop('check',fsetop('check',Data.(fnames{i})(im == j,:))');
  end
end
%-----------------------------------------------------------------------
function ok = l_comparehash(hash0,hash)
%L_COMPAREHASH Compares two hash structures.

ok = 1;
if hash0.nmonths ~= hash.nmonths
  warning('loadData:hash','Hash: Number of months mismatch.');
  ok = 0;
end
nmonths = max(hash0.nmonths,hash.nmonths);
if hash0.regions ~= hash.regions
  warning('loadData:hash','Hash: Region names or ordering mismatch.');
  ok = 0;
end

fnames = fieldnames(hash0);
if numel(fsetop('setxor',fnames,fieldnames(hash))) ~= 0
  warning('loadData:hash','Field names mismatch.');
  ok = 0;
end

for i = 1:numel(fnames)
  if any(strcmp(fnames{i},{'regions' 'nmonths'}))
    continue;
  end
  % gives a smooth behavior when nmonths increases:
  hsh0 = [hash0.(fnames{i}) zeros(1,nmonths-hash0.nmonths,'uint32')];
  hsh = [hash.(fnames{i}) zeros(1,nmonths-hash.nmonths,'uint32')];
  ind = find(hsh0 ~= hsh); % (now has the same length)
  if ~isempty(ind)
    warning('loadData:hash',sprintf(['Hash: field ''%s'' mismatch in index: [%s]'], ...
                    fnames{i},num2str(ind)));
    ok = 0;
  end
end

%-----------------------------------------------------------------------
function l_disp(quiet,str)
%L_DISP Display which can be turned on/off.

if ~quiet, disp(str); end

%-----------------------------------------------------------------------
