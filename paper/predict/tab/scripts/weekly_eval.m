%% Evaluate the "past predictions" table from .tex


% save .tex tables
if ~exist('savetofile','var')
    savetofile=false;
end


[~,~,~,N,date_tab,TAB] = readWeeklyTable();
% if ~date_tab(1) == DATES(end-7)
%   warning('Have you updated the table "tables/tablePredEval.tex"?');
% end

TAB = cat(1,TAB(:,:,1), TAB(:,:,2));
TAB = 1 - TAB; % easier to read 68% vs ~68 instead of 68% ~ 0.24
caption = ['Frequency (\%%) of all weekly reported predictions' ...
           ' that fell inside of (68/95\%% CrI, 7 days ahead), evaluated on' ...
           ' the following week (N~=~' num2str(numel(N)) ').'];

% caption = ['Andel (\%%) av samtliga veckorapporters prediktioner' ...
%             ' som fallit utanför (68/95\%% KI, 7 dagars), evaluerat när' ...
%             ' data kommit in (N~=~' num2str(numel(N)) ').'];
label = ['tab:PastWeekly'];
colspec = ['l' 'r' 'r' 'r'];
%
% rowlabel68 = regionList;
% rowlabel68{1} = [rowlabel68{1} ' (68\%%)'];
% rowlabel95 = regionList;
% rowlabel95{1} = [rowlabel95{1} ' (95\%%)'];
% rowlabel = [rowlabel68 rowlabel95];

tableWeekly = arr2latex(TAB*100,{'3.0f'}, ...
  'collabel',{'' 'Hospital (H)' 'Intensive (W)' 'Death (D)'},...
  'rowlabel',{'Uppsala (68\%%)' 'Stockholm' 'Sweden' ...
              'Uppsala (95\%%)' 'Stockholm' 'Sweden'}, ...
  'colspec',colspec, ...
  'hline','off','centering','on', ...
  'caption',caption, ...
  'label',label,...
  'h2line',3);


savepath = mfilename('fullpath');
savepath = savepath(1:end-11);

tabname = [savepath '../tableWeekly.tex'];
if savetofile
  fileID = fopen(tabname,'w');
  fprintf(fileID,'%s\n',tableWeekly);
  fclose(fileID);
  disp(['saved table: ' tabname]);
else
  disp(['didn''t save table: ' tabname]);
end