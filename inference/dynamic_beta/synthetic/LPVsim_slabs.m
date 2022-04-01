function x = LPVsim_slabs(F,x0,parstruct,slab)
%LPVSIM_SLABS Simulate linear parameter varying system
%   x = LPVSIM_SLABS(F,x0,beta,slab) simulates the system x_{k+1}=F_k*x_k,
%   x_0=x0.

% H. Runvik 2021-04-21

nstate=size(F{1},1);
ntime=size(parstruct.vals,1);
x=zeros(nstate,ntime);
xx=x0;

slabloop = 1+[0 find(diff(slab) ~= 0) numel(slab)];
slabloop = [slabloop(slabloop<=ntime) ntime+1];
nslab = numel(slabloop)-1;
for m = 1:nslab
    Fslab = F{m+slab(1)-1};
    for n=slabloop(m):slabloop(m+1)-1
        Fslab(parstruct.idx)=parstruct.vals(n,:);
        xx=Fslab*xx;
        x(:,n)=xx;
    end
end

end

