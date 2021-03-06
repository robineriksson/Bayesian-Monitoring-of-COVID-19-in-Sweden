% ERRORBARSREGION is a script that visualizes the similarity or
% discrepencies between the marginal posterior found for the different
% regions. The posterior in the figure is the mean and the standard
% deviation.
%
% The selection and representation of the marginals are the
% same as in PRIOR_POSTERIOR.m. If the script is called two times in a
% row, then the posterior values are retrieved from memory. Forcing
% new simulations can be done by 'clear ratenames'.

% Robin Marin (Eriksson) 2020-03-07
if ~exist('savetofile','var')
    savetofile=false;
end

if ~exist('verb','var')
    verb=false;
end


%clear ratenames
% load data
if ~exist('ratenames','var')
    if verb
        disp('reading from file');
    end
    ratenames = {'sigma', 'gammaI', 'gammaH' 'gammaW', ...
                 'R0' 'thetaA_' 'thetaE_','half_life',...
                 'E2I','HOSP','IC_HOSP','IFR'  };
    inverted = {'sigma', 'gammaI', 'gammaH' 'gammaW'};

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
    ending = '_1_100';
    [meanpost, stdpost, regionList] = l_getpost(postdate,ratenames,ending,inverted);

else
    if verb
        disp('using in-memory');
    end
end

% bar plot
regionid = 1:numel(regionList);
figure(1), clf

meanpost_plot = meanpost;
stdpost_plot = stdpost;

blue = [0 0 1]/256;
red = [1 0 0]/256;
for k = 1:12
    subplot(3,4,k)

    if ismember(ratenames{k}, {'E2I','HOSP','IC_HOSP','IFR'})
        meanpost_plot(k,:) = meanpost(k,:)*100;
        stdpost_plot(k,:) = stdpost(k,:)*100;
    else
        meanpost_plot(k,:) = meanpost(k,:);
        stdpost_plot(k,:) = stdpost(k,:);
    end



    bar(regionid, meanpost_plot(k,:));
    hold on
    er = errorbar(regionid,meanpost_plot(k,:), 1*stdpost_plot(k,:), -1*stdpost_plot(k,:));
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    title(titles{k},'interpreter','latex');
    xticks(regionid);

    ylim([max(0,min(meanpost_plot(k,:)-3*stdpost_plot(k,:))) ...
          max(0,max(meanpost_plot(k,:)+3*stdpost_plot(k,:)))])



    xticklabels(regionList);
    xtickangle(90);


    ax = gca;
    ax.XAxis.FontSize = 5;
    ax.TickLabelInterpreter = 'latex';

end

% polish target output size
set(gcf,'PaperPositionMode','auto');
%%
set(gcf,'Position',[100 100 850 500]);
%%
printpath = mfilename('fullpath');
printpath=[printpath(1:end-23) 'errorbars.pdf'];
if savetofile

    print('-dpdf',printpath);
    if verb
        disp(['saved figure: ' printpath]);
    end
else
    if verb
        disp(['Didn''t save figure: ' printpath]);
    end
end


function [meanpost, stdpost, regionList_tex] = l_getpost(postdate, ratenames, ending, inverted)
    %L_GETPOST(POSTDATE,RATENAMES,ENDING) get the posterior and return the mean and
    %   standard deviation for a specific date (POSTDATE) and rates
    %   (RATENAMES).
    %   ENDING   - {_1, _8, _15, _22}.
    %   INVERTED - name of rate that should be presented with its inverse.

    regionList = regions(false);

    regionList_tex = {'Stockholm' 'Uppsala' 'S\"{o}dermanland' '\"{O}sterg\"{o}tland' ...
                      'J\"{o}nk\"{o}ping' 'Kronoberg' 'Kalmar' 'Gotland' 'Blekinge' ...
                      'Sk\aa{}ne' 'Halland' 'V\"{a}stra G\"{o}taland' 'V\"{a}rmland' '\"{O}rebro' ...
                      'V\"{a}stmanland' 'Dalarna' 'G\"{a}vleborg' 'V\"{a}sternorrland' ...
                      'J\"{a}mtland' 'V\"{a}sterbotten' 'Norrbotten'}';

    file = [postpath() 'KLAM/perRegion/'];
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


function tag = l_biascomp(mu,sig,bias,ratenames)
    %L_BIASCOMP estimates if the bias is to large or not for the uncertainty
    %  given by 2*sig or not.

    bias_ = zeros(numel(ratenames),21);
    for i = 1:numel(ratenames)
        name = ratenames{i};
        rows = ismember(bias(:,2,1),name);
        if ismember(name,{'R0','IFR'})
            bias_(i,:) = mean(cell2mat(bias(rows,1,:)),1);
        else
            bias_(i,:) = cell2mat(bias(rows,1,:));
        end
    end


    for reg = 1:21
        d = 2*sig;
        tag = sqrt(bias_) >= 0.5*d;
    end
end