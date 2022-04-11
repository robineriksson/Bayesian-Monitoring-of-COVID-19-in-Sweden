%EVALSMOOTHING evaluates the smoothing and cleaning done by computing a
% by this scrupt defined metric. The evaluation is done per region on the
% C19 dataset.

datadate = [200401 210531] ;
data0 = loadData('C19');
tspan = find(data0.date == datadate(1)):find(data0.date == datadate(2));
data0.H = data0.H(tspan,:);
data0.W = data0.W(tspan,:);
data0.D = data0.D(tspan,:);
data0.Dinc = data0.Dinc(tspan,:);
data0.date=data0.date(tspan);

data1 = polishData(data0,'D','Dinc',1);
data1 = smoothData(data1,{'D' 'H' 'W'},{'Dinc' [] []});

HWD_1 = zeros(3,numel(tspan),21);
HWD_2  = HWD_1;
HWD_1(1,:,:) = data0.H;
HWD_1(2,:,:) = data0.W;
HWD_1(3,:,:) = data0.D;

HWD_2(1,:,:) = data1.H;
HWD_2(2,:,:) = data1.W;
HWD_2(3,:,:) = data1.D;

matr = squeeze(mean(max(abs(HWD_1-HWD_2)./max(1,abs(HWD_1)),[],1)));

eval_c = sprintfc('%g%%',round(100*matr,2,'significant'));
reg = regions();

meanMaxRelError = cat(2,reg,eval_c);

load Ncounties
weights = sum(N,1)./sum(N,[1 2]);
popweightmean = sprintfc('%g%%',round(100*weights*matr,2,'significant'));