
% nvar: max number of list
% fsort: field to sort

function ListVars(fpath,nvar,fsort)

if (nargin < 3),  fsort = '';  end
if (nargin < 2),  nvar = inf;  end


    vars = whos('-file',fpath);
    if strcmp(fsort,'bytes')
        [~,is] = sort([vars.bytes],'descend');
    else
        is = 1:length(vars);
    end

    for ii=1:min(nvar,length(is))
        fprintf('%s \t %10.2f MB \n',vars(is(ii)).name,vars(is(ii)).bytes/1e6); 
    end
    