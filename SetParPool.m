function SetParPool(nworker,minIdle)

if nargin < 2
    minIdle = 30;
end

p = gcp('nocreate');
if isempty(p)
    parpool(nworker,'IdleTimeout',minIdle);
else
    if p.NumWorkers ~= nworker
        delete(p);
        parpool(nworker,'IdleTimeout',minIdle);
    end
end
