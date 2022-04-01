%BOOTSTRAPTABLE checks the difference between two posterior files
% returns a table TAB which compares a posterior file produced from actual
% data with (x3) a synthetic set (generated using the Bayesian mean). For
% the synthetic set, we know the true parameters, and therefore the bias.
% Allowing for an imputed bias in for the actual data set.
%
% Columns in table format: [Region, bias^2, CoV, MSE, NMRSE]
%                                        
%
% R.Eriksson 2022-01-14

% save .tex tables
savetofile=false;



prefix = postpath;
prefix = [prefix 'SLAM/perRegion/'];
prefixSim = prefix;
ending = '_100';
if strcmp(ending,''), prefixSim = [prefix 'URDME/']; end
illustrate=false;
TAB = {'Region', 'bias$^2$','CoV [\%%]','NRMSE'};
marginal = {};
regions_all=regions();
for reg = 1:21
  region = regions_all{reg};
  fileData = [prefix 'slam210531_' region ...
    '_monthly_1' ending '.mat'];
  fileSim = {[prefixSim 'slam210531_' region ...
                     '_monthly_1_URDME1' ending '.mat'] ...
             [prefixSim 'slam210531_' region ...
                     '_monthly_1_URDME2' ending '.mat'] ...
             [prefixSim 'slam210531_' region ...
                     '_monthly_1_URDME3' ending '.mat']};
  filetrue = [prefix '../slam210531_mean_monthly_1'];
  [TAB_raw, marginal_raw] = posterror(fileData,fileSim,filetrue,illustrate);
  row = [region, TAB_raw(2,3), TAB_raw(3:end,2)'];
  TAB = cat(1,TAB,row);
  marginal = cat(3,marginal,marginal_raw);
end

caption = ['The Estimated bias$^2$ with Coefficient of Variation (CoV) '...
           ' and NRMSE, averaged over all parameters per region.'];
         
label = ['tab:bootstrap'];
colspec = ['l' 'r' 'r' 'r' 'r'];

% entries are actually string. We want them as numbers.
X = cell2mat(cellfun(@str2num,TAB(2:end,2:end),'un',0).')';
X = round(X,2,'significant'); % 2 significant points


% % add population per regional
% load Ncounties;
% Weights = sum(N(:,1:21))'/sum(N, [1 2]);
% X = cat(2,X,Weights);

% prepare latex output
tableBoot = arr2latex(X,{'g'}, ...
  'collabel',[TAB(1,:)],...% {'Population Weight'}],...
  'rowlabel',TAB(2:end,1),...
  'colspec',colspec, ...
  'hline','off','centering','on', ...
  'caption',caption, ...
  'label',label);




savepath = mfilename('fullpath');
savepath = savepath(1:end-14);
 
tabname = [savepath 'tableBootstrap.tex'];
if savetofile
  fileID = fopen(tabname,'w');
  fprintf(fileID,'%s\n',tableBoot);
  fclose(fileID);
  disp(['saved table: ' tabname]);
  
  marginal_name = [postpath() 'SLAM/marginal_bias.mat'];
  save(marginal_name,'marginal')
  disp(['saved marginals: ' marginal_name]);
else
  disp(['didn''t save table: ' tabname]);
end
