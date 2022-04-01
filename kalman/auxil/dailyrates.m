function rates = dailyrates(slabrates,slab)
%DAILYRATES Convert rates from per slab to per day
%   rates = DAILYRATES(slabrates,slab) returns a structure RATES with
%   identical fields as the input struct SLABRATES, except with a dalily
%   time resolution rather than slab-wise. The input SLAB defines the
%   length of the slabs.
%   
%  H. Runvik 2021-04-08

Nsamples = size(slabrates.sigma,2);
nslab = max(slab);
ntime=numel(slab);
fn = fieldnames(slabrates);
for k=1:numel(fn)
    f = slabrates.(fn{k});
    if isnumeric(f) && size(f,2)==Nsamples
        if size(f,1)==1
            rates.(fn{k})=repelem(f,ntime,1);
        elseif size(f,1)==nslab
            rates.(fn{k})=f(slab,:);
        else
            error('Incompatible slab size')
        end
    else
        rates.(fn{k}) = f;
    end
end

end

