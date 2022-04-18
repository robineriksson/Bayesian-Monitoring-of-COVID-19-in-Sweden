function [Uppsala,Stockholm,Sverige,week,date,TAB] = readWeeklyTable
%readWeeklyTable Parses the table tables/tablePredEval.tex.
%   [Uppsala,Stockholm,Sverige,week,date,TAB] = readWeeklyTable
%   returns the predictions for the 3 rows Hospital/Intensive
%   care/Death as (95%/68%/mean) confidence intervals in the first 5
%   columns [95% L, 68% L, mean, 68% U, 95% U]. The last column is the
%   actual outcome.
%
%   The vector week/date contains the number of the week/the date.
%
%   TAB is a frequency table for prediction misses. TAB is 3-by-3-by-3
%   for {Uppsala,Stockholm,Sverige}-by-{H,W,D}-by-68%/95 CI and is the
%   frequency of prediction misses.

% S. Engblom & H. Runvik 2021-05-16

% allocate output
week = [];
date = [];
Uppsala = {};
Stockholm = {};
Sverige = {};

% open file
savepath = mfilename('fullpath');
savepath = savepath(1:end-15);
fid = fopen([savepath 'tables/tablePredEval.tex']);

% read lines until first week was found
tline = fgetl(fid);
while numel(week) == 0
  if contains(tline,'\item[v.')
    k = strfind(tline,'v.');
    week = [week str2double(tline(k+2:end-1))];
  end
  tline = fgetl(fid);
end

% now process the rest
while ischar(tline)
  if contains(tline,'\item[v.')
    k = strfind(tline,'v.');
    week = [week str2double(tline(k+2:end-1))];
  end
  if contains(tline,'utfall')
    k = strfind(tline,'&');
    l = strfind(tline,'2');
    date = [date str2double(tline(l(1):k(1)-1))];
  end
  if contains(tline,'Uppsala') && contains(tline,'\CI')
    Uppsala{end+1} = l_textableread(tline);
  end
  if contains(tline,'Stockholm') && contains(tline,'\CI')
    Stockholm{end+1} = l_textableread(tline);
  end
  if contains(tline,'Sverige') && contains(tline,'\CI')
    Sverige{end+1} = l_textableread(tline);
  end
  tline = fgetl(fid);
end
fclose(fid);

% determine prediction accuracy
TAB = zeros(3,3,2); % {Uppsala,Stockholm,Sverige}-by-{H,W,D}-by-68%/95 CI
for ii = 1:numel(week)
  i = 0;
  for T = [Uppsala(ii) Stockholm(ii) Sverige(ii)]
    i = i+1;
    T_ = T{1};
    for j = 1:3 % H, W, D
      if prod(T_(j,[2 4])-T_(j,end)) > 0 % within 68%?
        TAB(i,j,1) = TAB(i,j,1)+1; % no, increment
        if prod(T_(j,[1 5])-T_(j,end)) > 0 % within 95%?
          TAB(i,j,2) = TAB(i,j,2)+1; % no, increment
        end
      end
    end
  end
end

% normalize to frequency
TAB = TAB/numel(week);
% (note: permute(mean(TAB,2),[1 3 2]) gives the accuracy regardless of {H,W,D})

% ----------------------------------------------------------------------
function data = l_textableread(tline)
%L_TEXTTABLEREAD Reads a single line of confidence intervals.

k1 = strfind(tline,'\CI');
k2 = strfind(tline,'(');
k3 = strfind(tline,') ');
data = zeros(3,6);
for i = 1:numel(k1)
  dataline = tline(k1(i)+4:k2(i)-4); % num}{num}{num}{num}{num
  newStr = split(dataline,'}{'); % num num num num num
  data(i,1:5) = cellfun(@str2num,newStr);
  data(i,6) = str2double(tline(k2(i)+1:k3(i)-1));
end
% ----------------------------------------------------------------------
