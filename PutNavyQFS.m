		if waitforbuttonpress
			cc = get(gcf,'CurrentCharacter');
			if strcmp(cc,'q')
				break;
			elseif strcmp(cc,'f')
				id = min(id+1,nd);
			elseif strcmp(cc,'s')
				id = max(id-1,1);
			end
		end
