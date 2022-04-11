%IFRtable computes the (bi)monthly IFR from the posterior

% R. Marin   2022-03-23 (Dinc weighting)
% R. Erisson 2021-10-19 (adjusted format)
% S. Engblom 2021-05-13
savetofile = false;
bimonthly = true; % per 2 months instead of per one

% Compare with URDME runs to estimate the robustness of the CrI
TOL_iqr = 0.32; % tolerances
URDME_k = [1];% [1 2 3], repeat with multiple bootstrap realizations

% Other than the (bi)monthly estimates, we want to find averages over
% these full intervals
intervals = {[200401 200501],... Initial month
             [200401 201101],... Up until start of November
             [200501 201101],...
             [200601 201101],... Up until start of November
             [200401 210531]}; % All data

% Weight the slabs (by Dinc) when averaging multiple ones into one estimate.
Dweight = false;
% Weight the regions [for the interval averages] on the size of the regions
Rweight = true;
% exclude transients
excludetransient = true;

%% Start the calculations
rng(0) % reproducability
abspath = mfilename('fullpath');
prefix = ['SLAM/perRegion/'];

regionList = regions;

postdate = '210531';

if Dweight % only load data if needed for weighting
  Data = loadData('C19');
  Data = polishData(Data); % removes negative deaths
  tspan_post = find((Data.date == 200401)):find((Data.date == str2num(postdate)));
  Dinc = Data.Dinc(tspan_post,:); % use the same data as in the posteior
end

% posterior filename
regname = [regionList(1:21); 'Sweden'];
ending_ = '_100'; % fast: '', slow full: '_100';
if strcmp(ending_, '_100')
  prefix_URDME = '';
else
  prefix_URDME = 'URDME/';
end
ending = ['_1' ending_];
post = strcat('slam',postdate,'_',regname,'_monthly', ending, '.mat');
filename = strcat(postpath,prefix,post);

% Compute the posterior weights when aggregating
Nc = load('Ncounties');
reglist = [1 2 22];

% for each URDME comparison we want to repeat the below code.
TAB_outer = cell(0,0);
for k = URDME_k
  ending_URDME = ['_1_URDME' num2str(k) ending_];

  post_URDME = strcat(prefix_URDME,'slam',postdate,'_',regname,'_monthly', ending_URDME, '.mat');
  filename_URDME = strcat(postpath,prefix,post_URDME);


  % seperate the replicate output and the true posterior
  TAB = cell(0,0);
  TAB_boot = cell(0,0);
  intervals_quantiles = zeros(0,1+1+5+2); %reg + mean + quantiles + dates
  for reg = reglist


    if reg == 22
      Weights = sum(Nc.N,1);
      file_URDME = filename_URDME(1:21);
      file = filename(1:21);
    else
      Weights = 1;
      file_URDME = filename_URDME(reg);
      file = filename(reg);
    end
    P_URDME = posteriorenger([],file_URDME,Weights,true);
    P = posteriorenger([],file,Weights,true);

    if Rweight
      Weights_R = Weights;
    else
      Weights_R = ones(size(Weights));
    end
    P_ = posteriorenger(1e5,file,Weights_R,false); % do a weighted fixed sample size draw

    % Perform below once for the URDME replicate and once on the actual
    % posterior
    penger = {P_URDME, P};
    for kk = [1 2]
      post = penger{kk};
      % assign each date a month
      [months,~,ixu] = fsetop('unique',...
        floor(post.meta.date/100)');

      % assign each month a mix of slabs (num months-by-num slabs)
      C = fsparse(ixu,post.meta.slabs,1);
      C = C./sum(C,2);
      % (each row gives the blend of slabs for the corresponding month)

      % need to remove first month (and likely also the last)
      if bimonthly
        % *** might need to change here to get the desired paring!
        remove = [1:5];
        months(remove) = []; % (for example, keep Mach)
        C(remove,:) = [];
        % if odd, remove another one from the beginning
        if mod(size(C,1),2)
          months(1) = [];
          C(1,:) = [];
        end
      end




      Nsamples = size(post.IFR,2);

      IFR = zeros(size(C,1),5); % IFR per month in the form of a 68%/95% CI
      IFRsample = nan(size(C,1),Nsamples); % (NaN is a safety measure)
      for i = 1:size(C,1)
        % sample from the slabs that goes into the mix
        [~,slabs,N] = find(C(i,:));
        % determine # samples per such slab
        N = round(N*Nsamples);
        % sample accordingly
        IFRsample_ = zeros(1,0);
        for j = 1:numel(slabs)
          % if not even sampling between months; when slabs span multiple
          % months.
          if Nsamples < N
            iN = randi(Nsamples,1,N(j));
          else
            iN = 1:N(j);
          end
          IFRsample_ = [IFRsample_ post.IFR(slabs(j),iN)];
        end
        IFRsample(i,:) = IFRsample_;
      end

      % pack up two and two if required
      if bimonthly
        IFRsample = [IFRsample(1:2:end,:) IFRsample(2:2:end,:)];
        months = months(1:2:end); % select 1st month in the pair
      end

      % pick out the stats median (NB!) and 68%/95% CI
      IFR = quantile(IFRsample',[0.05/2 0.32/2 0.5 1-0.32/2 1-0.05/2])';



      % last few months
      % im = find(2005 <= months & months < 2105);
      im = find(2004 <= months & months < 2106);

      %
      TAB_ = sprintf('\\CI{%4.2f}{%4.2f}{%4.2f}{%4.2f}{%4.2f}\n',...
        round(IFR(im,:),2,'significant')'*1e2);
      TAB_(end) = '';
      TAB_ = split(TAB_,newline);
      if kk == 1 % store
        TAB_boot = [TAB_boot TAB_];
      elseif kk == 2
        TAB = [TAB TAB_];
      end

    end

    % collect the average for the interval dates
    for j = 1:numel(intervals)
      % first find what slabs we want to look at.
      sel = find(ismember(P.meta.date(P.meta.slabstop),intervals{j}));
      sel = sel-1; % slabstop always includes 1, remove it.

      sel = sel(1):sel(2);
      if sel(1) == 0, sel=1;end
      val = P_.IFR(sel,:);



      % if Dinc weighted slabs
      if Dweight
        % Count the incidence death per slab
        if reg < 22
          Dinc_reg = Dinc(:,reg);
        else
          Dinc_reg = sum(Dinc,2);
        end

        % Find the cumulative Dinc per slab
        Dinc_cum = zeros(1,numel(P.meta.slabstop)-1);
        for i = 2:numel(P.meta.slabstop)
          Dinc_cum(i-1) = sum(Dinc(P.meta.slabstop(i-1):P.meta.slabstop(i)));
        end

        % Give zero weight to excluded slabs, and partial weight to partial slabs.
        Dinc_cum = Dinc_cum(sel) / sum(Dinc_cum(sel));
        val = Dinc_cum*val; % Weight the slabbed IFR
      end


      % quantile
      quant = [mean(val(:)) quantile(val(:), [0.05/2 0.32/2 0.5 1-0.32/2 1-0.05/2])];
      % output: [region id, mean, quantiles, dates]
      intervals_quantiles = cat(1, intervals_quantiles, [reg quant*100, intervals{j}]);
     end
  end


  % name of month + region
  monthName = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' ...
    'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
  if bimonthly
    monthName = strcat(monthName(:)','+',monthName([2:end 1]));
  end

  TAB = [monthName(mod(months(im),100))' TAB];
  TAB = [{''} regname(reglist)'; TAB];

  TAB_outer = [TAB_outer;TAB];


  % Now add notations
  % (asterisk): if |m_true - m_boot|/m_true > TOL_bias
  ifr_t = l_unwrap(TAB(2:end,2:end));
  ifr_b = l_unwrap(TAB_boot);
  asterisk = l_IQRcomp(ifr_t,ifr_b,[2 4], TOL_iqr);
  % (dagger): if bias >= 0.5*diam(A)
  dagger = l_biascomp(ifr_t,ifr_b,[2 4]);

  TAB_tex = l_texify(TAB_outer,asterisk, dagger, TOL_iqr)
end
%%
interval_averages = l_intervalwrap(intervals_quantiles, regname, Dweight, Rweight)
%%
%%
% save .tex tables
tabname = [abspath(1:end-13) '/../tableIFRBootstrap.tex'];
if savetofile
  %fileID = fopen(tabname,'w');
  [fileID, message] = fopen(tabname,'w');
  if fileID < 0
   error('Failed to open myfile because: %s', message);
  end
  fprintf(fileID,'%s\n',TAB_tex);
  fclose(fileID);
  disp(['saved table: ' tabname]);
else
  disp(['didn''t save table: ' tabname]);
end
return;
%%
figure(1), clf, hold on
for k = 1:size(P_.IFR,1)
  ifr = P_.IFR(k,:)*100;

  ifr_q = quantile(ifr,[0.5 0.025 0.975]);
  iqr = diff(ifr_q(2:3));
  disp([num2str(k) ': ' num2str(ifr_q) ', iqr: ' num2str(iqr)])
  [y,x] = ksdensity(ifr);
  plot(x,y)
end
legend(cellfun(@(x)num2str(x),num2cell(1:size(P_.IFR,1)),'UniformOutput',false))
hold off;

%%
figure(1), clf, hold on

  for j = 1:numel(intervals)
      % first find what slabs we want to look at.
      sel = find(ismember(P.meta.date(P.meta.slabstop),intervals{j}));
      sel = sel-1; % slabstop always includes 1, remove it.
      if sel(1) == 0, sel(1)=1;end
      sel = sel(1):sel(2);
      %if sel(2) == 0, sel=1;end
      val = P_.IFR(sel,:);

      [y,x] = ksdensity(val(:));
      plot(x,y)
  end
  legend(cellfun(@(x)num2str(x),intervals,'UniformOutput',false))
hold off

%% functions ...
function tex = l_texify(TAB, asterix, dagger, TOL_iqr)
%L_TEXIFY adds necessary additions into table TAB making it a bare minimum
%   table in TeX.
% ASTERIX are notes added on cells if certain conditions are
% meet. TOL_IQR is the defined tolerances used to generate
% the ASTERIX matrix.
%

tex = '';
tex = [tex '\begin{table}[htp]' newline ...
  '\centering' newline ...
  '\caption{' ...
  'Bi-monthly estimated IFR [\%] with 68\% CrI from October' newline ...
  '2020 to May 2021.}' newline ...
  '\label{tab:IFR}' newline ...
  '\begin{tabular}{l ' repmat('r ',1,size(TAB,2)-1) '}' newline ...
  '\hline'];
% per row
for row = 1:size(TAB,1)
  ROW = '';
  % per column, append a '&' sign
  for col = 1:size(TAB,2)
    val = TAB{row,col};

    if row > 1 & col > 1 & asterix(row-1,col-1)
      val = ['\footnotemark[3]' val];
    end
    if row > 1 & col > 1 & dagger(row-1,col-1)
      val = ['\footnotemark[4]' val];
    end



    if col == size(TAB,2) % final column, add '\\' instead of '&'
      ROW = [ROW val ' \\'];
    else
      ROW = [ROW val '& '];
    end

  end % end cold
  tex = [tex newline ROW];
  if row == 1 % after first row, add additional line
    tex = [tex '\hline'];
  end

end % end row

tex = [tex newline ...
      '\hline' newline ...
      '\end{tabular}' newline ...
      '\vspace{-0.25cm} \begin{flushleft}' newline ...
      '\begin{tabular}{l}' newline ...
      '{\tiny \footnotemark[4] The estimated bias is large compared' newline ...
      'to the 68\% CrI (see \matmet).}\\' newline ...
      '{\tiny \footnotemark[3] The 68\% CrIs of the posterior and the bootstrap' newline ...
      [' replicate do not share at least a ' num2str(100*(1-TOL_iqr)) '\%'] newline ...
      'overlap (see}\\' newline ...
      '\hspace{0.2cm}{\tiny \matmet).}' newline ...
      '\end{tabular}' newline ...
      '\end{flushleft}' newline ...
      '\end{table}\vspace{-0.5cm}'];

end

function val = l_unwrap(TAB)
%L_UNWRAP removews the \CI wrapping and returns a matrix with the entries
[row,col] = size(TAB);
val = zeros(row,5,col);
for r = 1:row
  for c = 1:col
    start = 5;
    for k = 1:5
      stop=start+3;
      val(r,k,c) = str2double(TAB{r,c}(start:stop));
      start = stop+3;
    end

  end
end
end

function [tag, score] = l_IQRcomp(val1,val2,cols,tol)
%L_IQRCOMP computes a (difference) score when comparing two IQRs


  % overlap x AND y (at the same time)
  set_and = @(X,Y)[max([X(:,cols(1),:) Y(:,cols(1),:)],[],2), ...
                   min([X(:,cols(2),:) Y(:,cols(2),:)],[],2)];

  % metric on the length, 0 = empty interval
  diam = @(X) max([diff(X,1,2), zeros(size(X))],[],2);

  x = squeeze(diam(set_and(val1,val2)));
  y = squeeze(diam(val1(:,cols,:))+diam(val2(:,cols,:)));

  score = 1-x./(0.5*y);
  tag = score > tol;
end

function tag = l_biascomp(val1,val2,cols)
%L_BIASCOMP estimates if the bias is to large or not for the uncertainty
%  given by val1(:,cols) or not.
  % metric on the length, 0 = empty interval
  diam = @(X) max([diff(X,1,2), zeros(size(X))],[],2);
  d = squeeze(diam(val1(:,cols,:)));
  bias = abs(squeeze(val1(:,3,:) - val2(:,3,:)));
  tag = bias >= 0.5*d;
end

function out = l_intervalwrap(interval, regname, Dweight, Rweight)
%L_INTERVALWRAP beautifies the interval output for presentation
out = cell(size(interval,1)+1,1);
% column ordering: reg, mean, qlolo qlohi qm qhili qhihi, date(1), date(2)
out{1} = sprintf('Dinc weighted: %d, Region weighted: %d',Dweight, Rweight);
for k = 1:size(interval,1)
  out{k+1} = sprintf('%s [%d-%d] \t %0.4f%% [%0.4f, %0.4f] (95%% CrI)', ...
    regname{interval(k,1)}, ... % region name
    interval(k,end-1), interval(k,end),... % dates
    round(interval(k,2),2,'significant'), ... % mean
    round(interval(k,3),2,'significant'), round(interval(k,end-2),2,'significant')); % 95% CrI
end

end
