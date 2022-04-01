function data = smoothData(data,cumfield,incfield)
%DATA = SMOOTHDATA(DATA,CUMFIELD,INCFIELD) mitigate reporting lag.
%   DATA = SMOOTHDATA(DATA,CUMFIELD,INCFIELD) smoothes the cumulative
%   fields CUMFIELD and the incidense fields INCFIELD to mitigate the
%   effects of lag or misses of reporting.
%
%   The algorithm removes unlikely incidenses in INCFIELD and spread
%   them out during the period immediately before in a way to mitigate
%   weekday dependencies in reporting. CUMFIELD is updated
%   accordingly.
%
%   When INCFIELD is empty, a simpler smoothing is performed on
%   CUMFIELD directly. Here a shape preserving Hermite interpolation
%   to the next different measurement is made, i.e., identical
%   measurements in a row are understood as missing reports and
%   replaced with a smoother transition.
%
%   Currently, the smoothed output fields are not integers.
%
%   Example:
%     Data = loadData('C19');
%     Data = polishData(Data,'D','Dinc');
%     Data0 = Data;
%     Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});
%     reg = 1; % region to inspect...
%
%     % inspect smoothing of cumulative fields
%     figure(1), clf,
%     tspan = 1:size(Data.date);
%     stairs(tspan,Data0.H(:,reg),'g.-'); hold on
%     stairs(tspan,Data.H(:,reg),'b.-');
%     stairs(tspan,Data0.W(:,reg),'g.-','Handlevisibility','off');
%     stairs(tspan,Data.W(:,reg),'b.-','Handlevisibility','off');
%     legend('Original','Update');
%
%     % pull out data for inspection
%     figure(2), clf,
%     ixnotnan = find(~isnan(Data.Dinc(:,reg)));
%     date = Data.date(ixnotnan(1));
%     Dinc = Data.Dinc(ixnotnan,reg);
%     Dinc0 = Data0.Dinc(ixnotnan,reg);
%     tspan = (1:numel(Dinc))';
%
%     s = num2str(date);
%     s = [s(1:2) '-' s(3:4) '-' s(5:6)];
%     wd = mod(weekday(s)-1+(0:numel(tspan)-1),7)+1;
%
%     Dinc_smooth = max(sgolayfilt(Dinc,1,7),0);
%     Dcomp_inc = poissrnd(Dinc_smooth);
%
%     plot(tspan,Dinc0,'g', ...
%          tspan,Dinc,'b', ...
%          tspan,Dcomp_inc,'r--');
%     legend('Original','Update','Synthetic');
%
%     figure(3), clf,
%     d = Dinc-Dinc_smooth;
%     Z = abs(d)./max(sqrt(Dinc_smooth),1);
%
%     Dinc_smooth0 = max(sgolayfilt(Dinc0,1,7),0);
%     d0 = Dinc0-Dinc_smooth0;
%     Z0 = abs(d0)./max(sqrt(Dinc_smooth0),1);
%
%     Dcomp_inc_smooth = max(sgolayfilt(Dcomp_inc,1,7),0);
%     Zcomp = abs(Dcomp_inc-Dcomp_inc_smooth)./ ...
%        max(sqrt(Dcomp_inc_smooth),1);
%     semilogy(tspan,sort(Z0),'g', ...
%        tspan,sort(Z),'b',tspan,sort(Zcomp),'r--');
%     ylim([1e-1 1e2]);
%     legend('Original','Update','Synthetic');
%     xlabel('Normalized difference (Z)');
%
%     figure(4), clf,
%     bar([sparse(wd,1,Dinc0)./sparse(wd,1,1)/sum(Dinc0) ...
%          sparse(wd,1,Dinc)./sparse(wd,1,1)/sum(Dinc)  ...
%          sparse(wd,1,Dcomp_inc)./sparse(wd,1,1)/sum(Dcomp_inc)]);
%     set(gca,'xticklabel',{'Sun' 'Mon' 'Tue' 'Wed' 'Thu' 'Fri' 'Sat'});
%     legend('Original','Update','Synthetic');
%
%     figure(5), clf,
%     stairs(Data0.D(:,reg),'g'); hold on,
%     stairs(Data.D(:,reg),'b'); hold on,
%     legend('Original','Update');
%
%   See also POLISHDATA, LOADDATA.

% S. Engblom 2021-05-01

if nargin < 3
  incfield = {};
  if nargin < 2
    cumfield = {};
  end
end
if ~iscell(cumfield), cumfield = {cumfield}; end
if ~iscell(incfield), incfield = {incfield}; end
if numel(cumfield) ~= numel(incfield)
  error('Can not match cumulative and incremental fields.');
end

% add a note that data has been smoothed
if contains(data.rev,'smoothData')
  warning('Data seems to have already been smoothed once.');
end
data.rev = [data.rev ' | smoothData'];

% parameter Ndays: time window backwards within which the fix will
% (recusively) take place
Ndays = 28;
for ix = 1:numel(incfield)
  % if incfield is empty, smooth only cumfield
  if isempty(incfield{ix})
    for reg = 1:size(data.(cumfield{ix}),2)
      D = data.(cumfield{ix})(:,reg);
      ixnan = find(~isnan(D));
      D = D(ixnan);
      tspan = 1:numel(D);
      D_inc = diff(D);
      % shape preserving Hermite interpolation to next "suitable"
      % measurement:

      % "very smooth":
% $$$       id = find(D_inc ~= 0);

      % more careful than above: don't smooth small jumps and also turn off
      % smoothing whenever both 'before' and 'after' the jump are
      % small, additionally don't ramp up too slowly (length of ramp
      % and jump should be at least equal)
      id = find(D_inc ~= 0);
      jump = abs(D_inc(id));     % interpolation jumps
      % no interpolation for small jumps (also takes away jumps
      % between small states!)
      jump(jump <= 5) = 1;
      len = diff([id(1)-1; id]); % ditto lengths
      % "warning of jump" too early?
      missing = len > jump;
      % then add a point closer to the jump:
      missing = id(missing)-floor(jump(missing));
      id = unique([id; missing]);

      % finalize interpolation indices
      id = [(1:id(1))'; id+1];
      if id(end) ~= numel(D)
        id = [id; numel(D)]; % add end if missing
      end
      pp = pchip(tspan(id),D(id));
      Dsmooth = ppval(pp,tspan);
      data.(cumfield{ix})(ixnan,reg) = Dsmooth;
    end
  else
    % incidense-based smoothing
    for reg = 1:size(data.(incfield{ix}),2)
      % pull out (D,Dinc)
      date = data.date;
      D = data.(cumfield{ix})(:,reg);
      Dinc = data.(incfield{ix})(:,reg);
      ixnan = find(~isnan(Dinc));
      D = D(ixnan);
      Dinc = Dinc(ixnan);
      date = date(ixnan);

      % days of the week
      s = num2str(date(1));
      s = [s(1:2) '-' s(3:4) '-' s(5:6)];
      wd = mod(weekday(s)-1+(0:numel(date)-1),7)+1; % 1:7
      
      % "peel" layers of outliers and remove with smoothed value
      for i = 1:9 % (why 9? see below)
                  % construct smooth baseline for comparision purposes
        Dinc_smooth = max(sgolayfilt(Dinc,1,7),0);

        % find all points to within (11-i)*SD of the baseline, where 1SD =
        % sqrt(baseline) (Poissonian assumption)
        d = Dinc-Dinc_smooth;
        ixoutlier = find(abs(d) > (11-i)*max(sqrt(Dinc_smooth),1));
        % (hence the loop goes from +/- 10SD... +/- 2SD, i.e., 95% CI)

        % replace those points with the smooth value
        Dinc2 = Dinc;
        Dinc2(ixoutlier) = Dinc_smooth(ixoutlier);

        % distribute the indeuced diff backwards @ a fixed rate for a total of
        % Ndays days
        delta = d(ixoutlier);
        ixbwd = max(ixoutlier+(-1:-1:-Ndays),1);
        % rate model: rate per day during Ndays days
        wc = full(sparse(wd,1,Dinc)); % count per weekday
        wc = 7*wc/sum(wc); % normalized to unit rate per day
        % weekday count adjustment: "complement of measured" and, for negative
        % correction,"same as measured"
        rate = (delta >= 0).*max(2-wc(wd(ixbwd)),0)+ ...
               (delta < 0).*wc(wd(ixbwd)).*(Dinc2(ixbwd) > -delta);
        % (rates > 0, and = 0 when delta < 0 and Dinc2 seems small)
        rate = reshape(rate,size(ixbwd));
        % (depending on the sign of the diff to distribute)
        rate = rate.*(Ndays:-1:1); % backwards decay
        %rate = rate.*cumsum(1./(1:Ndays),2,'reverse');
        rate = rate./sum(rate,2); % normalization

        % distribute the delta accordingly
        Delta = delta.*rate;
        Dinc2 = Dinc2+sparse(ixbwd,1,Delta,numel(Dinc2),1);

        Dinc = Dinc2; % next round
      end

      % recreate D and Dinc
      D = cumsum(Dinc);
      data.(cumfield{ix})(ixnan,reg) = D;
      data.(incfield{ix})(ixnan,reg) = Dinc;
    end
  end
end
