function [init] = getInit(useinit, filepath, final)
%INIT = GETINIT(USEINIT, REGION, FINAL) Retrieves the initial guess for SLAM,
%   i.e., and old posterior and constructs the wanted structure of the
%   output. Similar to posteriorenger, however, with additional steps.
%
%   UNSEINIT - true or false, if nonempty or empty
%   FILEPATH - the filepath to the posterior to retrieve.
%   FINAL    - true or false, if true start from last value, otherwise
%              start from mean value.

% R. Eriksson 210522

if useinit
 init = struct();
 rates = posteriorenger([],filepath);

 % If old format is given, convert to new the new format.
 try
   rates.E2I = rates.F0;
   rmfield(rates,'F0');
 end
 try
   rates.A2I = rates.F1;
   rmfield(rates,'F1');
 end



 init.meta = rates.meta; % when made into mat, this is lost.
 [rates0, rates0_names] = struct2mat(rates);
 if final
  init.rates0 = rates0(:,end);
 else
  init.rates0 = mean(rates0,2);
 end
 init.ratesAll = rates0;
 try
   init.nslab = init.meta.slabs(end);
 catch
   init.nslab = init.meta.slab(end);
 end
 init.names = rates0_names;
else
 init=[];
end
end