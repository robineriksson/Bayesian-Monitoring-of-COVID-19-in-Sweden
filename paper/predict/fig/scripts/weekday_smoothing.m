%WEEKDAY_SMOOTHING illustrates the periodic effect that
% POLISHDATA and SMOOTHDATA has on the C19 data.

% R. Marin 2022-04-06 (formating of example in SMOOTHDATA (c) S.Engblom)

if ~exist('savetofile','var')
    savetofile = false;
end

if ~exist('reg','var')
    reg = 1;
end


Data = loadData('C19');

datetime = [200401 210531]
tspan_date = find(Data.date == datetime(1)):find(Data.date == datetime(2));
Data.date = Data.date(tspan_date);
Data.H = Data.H(tspan_date,:);
Data.W = Data.W(tspan_date,:);
Data.D = Data.D(tspan_date,:);
Data.Dinc = Data.Dinc(tspan_date,:);


Data = polishData(Data,'D','Dinc');
Data0 = Data;
Data = smoothData(Data,{'D' 'H' 'W'},{'Dinc' [] []});


tspan = 1:size(Data.date);



ixnotnan = find(~isnan(Data.Dinc(:,reg)));
date = Data.date(ixnotnan(1));
Dinc = Data.Dinc(ixnotnan,reg);
Dinc0 = Data0.Dinc(ixnotnan,reg);
tspan = (1:numel(Dinc))';


s = num2str(date);
s = [s(1:2) '-' s(3:4) '-' s(5:6)];
wd = mod(weekday(s)-1+(0:numel(tspan)-1),7)+1;

% URDME
load('URDME/URDMEoutput/URDME_all')
if D.date(1) ~= datetime(1) || D.date(end) ~= datetime(end)
    warning('start and/or stop dates are different in URDME')
end
D_URDME = squeeze(D.U(ixnotnan,7,reg,1))'; % use only 1 sample
Dinc_URDME = [0 diff(D_URDME)];
Dinc_smooth = max(sgolayfilt(Dinc_URDME,1,7),0);



figure(1), clf,
bar([sparse(wd,1,Dinc0)./sum(Dinc0) ...
     sparse(wd,1,Dinc)./sum(Dinc)  ...
     sparse(wd,1,Dinc_URDME)./sum(Dinc_URDME)].*100);
ylabel('$\%$','Interpreter','Latex');
set(gca,'xticklabel',{'Sun' 'Mon' 'Tue' 'Wed' 'Thu' 'Fri' 'Sat'})
ax = gca;
ax.TickLabelInterpreter = 'latex';

legend('Original','Smoothed','Synthetic','interpreter','latex');

% polish target output size
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 500 350]);

printpath = mfilename('fullpath');
printpath = [printpath(1:end-25) 'smoothing.pdf'];
if savetofile
    print('-dpdf',printpath);
    disp(['saved figure: ' printpath])
else
    disp(['didn''t save figure: ' printpath])
end