%% Navigate figures w/ waitforbuttonpress

% vfs: variable for FS
% mfs: [min max] of vfs
% ved: variable for ED
% med: [min max] of ved

% e.g.,             vfs = 'is';  mfs = [1 ns];  ved = '';  med = [];  NavigateFig;

        if waitforbuttonpress
			cc = get(gcf,'CurrentCharacter');
			if strcmp(cc,'q')
				eval('break');
			elseif length(vfs) > 0 && strcmp(cc,'f')
				eval([vfs ' = min(' vfs '+1,' num2str(mfs(2)) ');']);
			elseif length(vfs) > 0 && strcmp(cc,'s')
				eval([vfs ' = max(' vfs '-1,' num2str(mfs(1)) ');']);
			elseif length(ved) > 0 && strcmp(cc,'e')
				eval([ved ' = min(' ved '+1,' num2str(med(2)) ');']);
			elseif length(ved) > 0 && strcmp(cc,'d')
				eval([ved ' = max(' ved '-1,' num2str(med(1)) ');']);
			end
		end
