function printing(type, text, i, inum)
% this functions handles printing loop progress. 
	persistent FPRINTF_R_LENGTH;
	if type == 1
		str = sprintf('Finished %s calculations for %d of %d images.', text, i, inum);
		lstr = length(str);
		if(i==1)
			FPRINTF_R_LENGTH=fprintf(str);
		else
			if FPRINTF_R_LENGTH>0
				b_string = repmat(sprintf('\b'), 1, FPRINTF_R_LENGTH);
				fprintf(b_string);
			end
			FPRINTF_R_LENGTH=fprintf(str);
		end
	elseif type == 2
		fprintf('Finished %s calculations for %d of %d images.\r', text, i, inum);
	end

end