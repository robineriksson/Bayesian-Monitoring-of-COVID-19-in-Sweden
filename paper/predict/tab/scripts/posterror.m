function [TAB,marginal] = posterror(fileData,fileSim,fileTrue,illustrate)
%POSTERROR returns an error measurments, TAB, for a posterior with unknown
%   truth (FILEDATA) using a posterior with known truth (FILESIM) which has
%   the truth found in FILETRUE. ILLUSTRATE {false} also creates
%   comparative plots to visualize the difference between the two
%   posteriors.
%
%   The TAB includes:
%   CoV: Coefficient of variation, indicates the
%      spread of the distribution,
%   CoB: Coefficient of bias, swaps the
%      std for the bias in the CoV.
%   NRMSE: Normalized root mean square error. Mean square error, here
%      the induced bias from FILESIM is used for FILEDATA, resulting in a
%      induced error measurment.
%
%


% load posterior files
postData = posteriorenger([],fileData);
[postDataMean, postDataVar,dataX,names] = l_cellsum(postData);
region = postData.meta.region; % for plotting and if URDME is used.



for k = 1:numel(fileSim)
  try
    postSim = posteriorenger([],fileSim{k});
    [mmean, vvar,simX] = l_cellsum(postSim);
    if k == 1
      postSimMean = zeros(numel(mmean),numel(fileSim));
      postSimVar = zeros(numel(vvar),numel(fileSim));
      regionSim = postSim.meta.region; % for plotting and if URDME is used.
    end
    postSimMean(:,k) = mmean;
    postSimVar(:,k)  = vvar;

  catch
    warning(['did not find file' fileSim{k}]);
  end

end

if ~strcmp(region,regionSim)
  error('trying to compare two different regions');
end


% load "truth"
load(fileTrue,'rates');
postTrue = rates;
postTrue.R0 = squeeze(postTrue.R0(:,1,strcmp(regions,region)));
postTrue.beta = squeeze(postTrue.beta(:,1,strcmp(regions,region)));

[postTrueMean,~,trueX,fieldnames] = l_cellsum(postTrue); % remove not


% we want to infer the error in postData, for this we first compute the
% bias in postSim with postTrue (which should be similar, known truth).
bias2 = mean((postSimMean - postTrueMean).^2,2);
CoV   = 100*[sqrt(postDataVar) sqrt(mean(postSimVar,2))]./[postDataMean mean(postSimMean,2)];
CoB   = 100*[sqrt(bias2) sqrt(bias2)]./[postDataMean mean(postSimMean,2)];
MSE   = [postDataVar mean(postSimVar,2)] + bias2;
NRMSE = 100*sqrt(MSE) ./ [postDataMean mean(postSimMean,2)];

TAB = cell(0,2);
%meanfunc = @(x) mean(x);
%meanfunc = @(x) sqrt(mean(x.^2));
meanfunc = @(x) median(x);

TAB = cat(1,TAB,strsplit(num2str(meanfunc(CoV))));
TAB = cat(1,TAB,strsplit(num2str(meanfunc(CoB))));
TAB = cat(1,TAB,strsplit(num2str(meanfunc(NRMSE))));
TAB = [{'CoV [%]' 'CoB [%]' 'NRMSE [\%]' }' TAB]; % rownames
TAB = [{region 'Data' ['Simulated x ' num2str(size(postSimMean,2))] }; TAB];

marginal = cat(2,num2cell(bias2),fieldnames);

if illustrate
  figure(1),clf
  subplot(3,1,1)
  boxplot(postData.R0','Color','b','symbol', '')
  hold on;
  boxplot(postSim.R0','Color','r','symbol', '')

  plot(postTrue.R0,'k*')
  title(['R0 | ' region ', blue: data, '...
                        'red: simulated, '...
                        'black: truth (sim)']);

  hold off

  subplot(3,1,2)
  boxplot(postData.IFR','Color','b','symbol', '')
  hold on;
  boxplot(postSim.IFR','Color','r','symbol', '')
  plot(postTrue.IFR,'k*')
  title(['IFR | ' region ', blue: data, '...
                         'red: simulated, '...
                         'black: truth (sim)']);


  subplot(3,1,3)
  static = ~fsetop('ismember',names,{'R0' 'IFR'});
  boxplot(dataX(static,:)','Color','b')
  hold on
  boxplot(simX(static,:)','Color','r')
  plot(trueX(static,:),'k*')
  title(['static | ' region ', blue: data, '...
                            'red: simulated, '...
                            'black: truth (sim)']);


end

end

function [mu,sigma2,X,names] = l_cellsum(rates)
%L_CELLSUM computes the mean (MU) and the variance (SIGMA2) of the
%   different fields in RATES.

[X,names] = struct2mat(rates);


exclude = {'A2I' 'rho' 'thetaI' 'thetaI_' 'thetaA' 'thetaE'  ...
           'F0ave' 'F1ave' 'F2ave' 'F2dave' 'F3ave' 'F3dave' ...
           'F4ave' 'beta' 'meta' 'lambda' 'gammaA' ...
           'R0' 'IFR'};
include = ~fsetop('ismember',names,exclude);


names_ = names(include);

dynamic = {'R0' 'IFR'}';
dmu = []; dsigma2 = [];
for id = 1:numel(dynamic)
    dname = dynamic{id};
    dmu = [dmu; mean(X(ismember(names,dname),:),[1 2])];
    dsigma2 = [dsigma2; var(X(ismember(names,dname),:),[],[1 2])];
end
mu = mean(X(include,:),2);
sigma2 = var(X(include,:),[],2);


mu = [mu; dmu];
sigma2 = [sigma2; dsigma2];
names = cat(1,names_,dynamic);
end
