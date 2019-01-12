function SetParPool(nworker)

p = gcp('nocreate');
if isempty(p)
    parpool(nworker);
else
    if p.NumWorkers ~= nworker
        delete(p);
        parpool(nworker);
    end
end
