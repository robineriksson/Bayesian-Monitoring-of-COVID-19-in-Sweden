%POSTERIOR_CLUSTER is a script that visualizes the similarity or
%   discrepencies between the posterior found for different regions.

% Robin Marin (Eriksson) 2020-3-7
save2file = true;

%clear ratenames
% load data
if ~exist('ratenames','var')
  disp('reading from file');
  ratenames = {'sigma', 'gammaI', 'gammaH' 'gammaW', ...
               'R0' 'thetaA_' 'thetaE_','half_life',...
               'E2I','HOSP','IC_HOSP','IFR'  };
  inverted = {'sigma', 'gammaI', 'gammaH' 'gammaW'};
%   titles = {'$\sigma^{-1}$', '$\gamma_I^{-1}$', '$\gamma_H^{-1}$', '$\gamma_W^{-1}$', ...
%             '$R_t$' '$\theta_A$', '$\theta_E$','half life', ...
%             '$E_2I$','HOSP','IC HOSP', 'IFR'};
  titles= {'$\sigma^{-1}$ [days]' ...
                      '$\gamma_I^{-1}$ [days]' ...
                      '$\gamma_H^{-1}$ [days]'...
                      '$\gamma_W^{-1}$ [days]' ...
                      '$R_t$' ...
                      '$\theta_A$ [$A \rightarrow \varphi$]' ...
                      '$\theta_E$ [$E \rightarrow \varphi$]' ...
                      'Virus half-life [hours]' ...
                      '$E_2I$ [$E \rightarrow I$, \%]' ...
                      'HOSP [$I \rightarrow H$, \%]' ...
                      'IC HOSP [$H \rightarrow W$, \%]' ...
                      'IFR [$E \rightarrow D$, \%]'};
  postdate = '210531';
  ending = '_1';
  [meanpost, stdpost, regionList] = l_getpost(postdate,ratenames,ending,inverted);
  
else
  disp('using in-memory');
end

% bar plot
regionid = 1:numel(regionList);
tiledlayout(3,4)

for k = 1:12
  nexttile
%   if ~ismember(k,[5,12])
%     i = i+1;
%   
  %subplot(2,5,i)
  bar(regionid, meanpost(k,:));
  hold on
  er = errorbar(regionid,meanpost(k,:), 1*stdpost(k,:), -1*stdpost(k,:));
  er.Color = [0 0 0];
  er.LineStyle = 'none';
  title(titles{k},'interpreter','latex');
  xticks(regionid);
  
  ylim([max(0,min(meanpost(k,:)-3*stdpost(k,:))) ...
        max(0,max(meanpost(k,:)+3*stdpost(k,:)))])
        

      
  xticklabels(regionList);
  xtickangle(90);
  
  

  %set(gca,'fontsize',5)

  
 
  
  ax = gca;
  ax.XAxis.FontSize = 5;
  ax.TickLabelInterpreter = 'latex';
%   else
%     set(gca,'visible','off')
%   end
%     
end

% polish target output size
set(gcf,'PaperPositionMode','auto');
%%
set(gcf,'Position',[100 100 900 500]);
%%  
if save2file
  printpath = mfilename('fullpath');
  printpath=[printpath(1:end-23) 'errorbars.pdf'];
%   printpath = [printpath(1:end-10) '..' printpath(end-10:end)];
%   printpath = [printpath '_' strrep(region,' ','_')];
%   printpath = strrep(printpath,'å','a');
%   printpath = strrep(printpath,'ä','a');
%   printpath = strrep(printpath,'ö','o');
  %print('-depsc',printpath);
  print('-dpdf',printpath);
  disp(['saved figure: ' printpath]);
  % *** does not always work, run manually in this case:
  %unix(['epstopdf ' printpath '.eps']);
end


function [meanpost, stdpost, regionList_tex] = l_getpost(postdate, ratenames, ending, inverted) 
%L_GETPOST(POSTDATE,RATENAMES,ENDING) get the posterior and return the mean and
%   standard deviation for a specific date (POSTDATE) and rates
%   (RATENAMES).
%   ENDING   - {_1, _8, _15, _22}.
%   INVERTED - name of rate that should be presented with its inverse. 

regionList = {'Stockholm' 'Uppsala' 'Södermanland' 'Östergötland' ...
  'Jönköping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
  'Skåne' 'Halland' 'Västra Götaland' 'Värmland' 'Örebro' ...
  'Västmanland' 'Dalarna' 'Gävleborg' 'Västernorrland' ...
  'Jämtland' 'Västerbotten' 'Norrbotten'}';

regionList_tex = {'Stockholm' 'Uppsala' 'S\"{o}dermanland' '\"{O}sterg\"{o}tland' ...
  'J\"{o}nk\"{o}ping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
  'Sk\aa{}ne' 'Halland' 'V\"{a}stra G\"{o}taland' 'V\"{a}rmland' '\"{O}rebro' ...
  'V\"{a}stmanland' 'Dalarna' 'G\"{a}vleborg' 'V\"{a}sternorrland' ...
  'J\"{a}mtland' 'V\"{a}sterbotten' 'Norrbotten'}';

file = [postpath() 'SLAM/perRegion/'];
meanpost = zeros(numel(ratenames),numel(regionList));
stdpost = zeros(numel(ratenames),numel(regionList));

for i = 1:numel(regionList)
  region = regionList{i};
  try
    rates = posteriorenger([],...
      [file 'slam' postdate '_' region '_monthly' ending '.mat']);
    %load([file 'slam' postdate '_' region '_monthly' ending '.mat'],'rates');
    for j = 1:numel(ratenames)
      rate = rates.(ratenames{j});
      if any(strcmp(ratenames{j},inverted))
        rate = 1./rate;
      end
      meanpost(j,i) = mean(rate,[1 2]);
      stdpost(j,i) = std(rate,[],[1 2]);
    end
  catch
    warning([region ' was not found']);
    regionList_tex{i} = [];
  end
end
end

