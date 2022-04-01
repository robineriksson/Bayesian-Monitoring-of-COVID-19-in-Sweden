function [J, L] = betaCost_S(beta,y,F,H,betaIdx,wbeta,slab,S,x0)
%BETACOST_S cost function for ss model with time varying parameter beta

% H. Runvik 2021-01-12

nstate=size(F{1},1);
ntime=size(y,2);
if nargin<9
    x0=beta(1:nstate)';
    beta=beta(nstate+1:end);
end
L=0;
x=x0;

slabloop = 1+[0 find(diff(slab) ~= 0) numel(slab)];
slabloop = [slabloop(slabloop<=ntime) ntime+1];
nslab = numel(slabloop)-1;
for m = 1:nslab
    Fslab = F{m+slab(1)-1};
    for n=slabloop(m):slabloop(m+1)-1
        Fslab(betaIdx)=beta(n);
        x=Fslab*x;
        yy = y(:,n)-H*x;
        L = L -0.5*(yy'*(S(:,:,n)\yy));%+size(y,1)*log(2*pi)+log(det(S(:,:,n))));
    end
end

% F3=permute(reshape(cell2mat(F(slab))',nstate,nstate,[]),[2,1,3]);
% betaIdx3 = repmat(betaIdx,1,1,ntime);
% F3(betaIdx3) = beta;
% 
% F3cumprod =zeros(nstate,nstate,ntime);
% F3cumprod(:,:,1) = F3(:,:,1);
% for n=2:ntime % something faster possible?
%     F3cumprod(:,:,n) = F3(:,:,n)*F3cumprod(:,:,n-1);
% end
% allx = reshape(reshape(permute(F3cumprod,[2 1 3]),nstate,[]).'*x0,nstate,[]);
% ally = H*allx;
% diffy = y-ally;

% L=0;
% 
% for n=1:ntime % something faster possible?
%     L = L -0.5*(diffy(:,n)'*(S(:,:,n)\diffy(:,n)));
% end


J = 1e-4*(wbeta*sum(diff(beta).^2) - L);