function I = ipreproc(I, variables)
% this functions handles preprocessing for one image. 

    %%% More variables %%%
    meanIn = 127;
	H = 1;
	if isfield(variables, 'quality') == 1
		if(variables.quality.blur_size > 0)
		   H = fspecial('gaussian', variables.quality.blur_size, variables.quality.blur_sigma);
		end
	end

    %%% Get desired channel %%%
    if 1 == strcmp(variables.data.channel, 'gray')
        I = rgb2gray(I);
    elseif 1 == strcmp(variables.data.channel, 'red')
        I = I(:,:,1);
    elseif 1 == strcmp(variables.data.channel, 'green')
        I = I(:,:,2);
    elseif 1 == strcmp(variables.data.channel, 'blue')
        I = I(:,:,3);
    end

    % check if RGB image
    [channels] = size(I,3);
	
	%%% average the image over a window %%%
	if isfield(variables.data, 'mean_window') == 1
		if channels == 1 && variables.data.mean_window >= 3
			H = ones(variables.data.mean_window);
			I = imfilter(double(I),H,'replicate');
			I = uint8(I./(variables.data.mean_window*variables.data.mean_window));
		end
	end

    %%% resize %%%
	if ~strcmp(variables.patches.type,'park')
		I = imresize(I, [variables.data.max_y variables.data.max_x]);
	end
	
	%%% quality perterbation %%%
	if isfield(variables, 'quality') == 1
		if(variables.quality.blur_num > 0)
		   for q=1:variables.quality.blur_num
			  I = imfilter(I,H,'replicate');
		   end
		end
		%if(variables.quality.downsample < 1)
		%	I = imresize(I, [floor(variables.data.max_y*variables.downsample) floor(variables.data.max_x*variables.downsample)]); 
		%end
	end
	
	%%% histeq %%%
    if(variables.data.histeq == 1)
        if(channels~=3)
            I = histeq(I);
        else %histogram equalization to preserve color information
            I = double(I)/255;
            srgb2lab = makecform('srgb2lab');
            lab2srgb = makecform('lab2srgb');
            shadow_lab = applycform(I, srgb2lab); % convert to L*a*b*
            max_luminosity = 100;
            L = shadow_lab(:,:,1)/max_luminosity;
            shadow_histeq = shadow_lab;
            shadow_histeq(:,:,1) = histeq(L)*max_luminosity; %perform histeq on luminance channel
            I = applycform(shadow_histeq, lab2srgb); %transform back to rgb
            clearvars srgb2lab lab2srgb shadow_lab max_luminosity L shadow_histeq;
            %place code here to implement different colorspaces :) maybe use variables.channel to set this?
        end
    end

    %%% ellipse mask %%%
    if ~(variables.data.ellipse_a == 0 && variables.data.ellipse_b == 0)
        if (channels~=3)
			if (variables.data.ellipse_type == 1)
				I = ellipseMatrix(variables.data.ellipse_y, variables.data.ellipse_x, variables.data.ellipse_a, variables.data.ellipse_b, variables.data.ellipse_t, I, meanIn);
			elseif (variables.data.ellipse_type == 2)
				I = ellipseMatrixInv(variables.data.ellipse_y, variables.data.ellipse_x, variables.data.ellipse_a, variables.data.ellipse_b, variables.data.ellipse_t, I, meanIn);
			end
        else%place mask on all three channels
            I(:,:,1) = ellipseMatrix(variables.data.ellipse_y, variables.data.ellipse_x, variables.data.ellipse_a, variables.data.ellipse_b, variables.data.ellipse_t, I(:,:,1), meanIn);
            I(:,:,2) = ellipseMatrix(variables.data.ellipse_y, variables.data.ellipse_x, variables.data.ellipse_a, variables.data.ellipse_b, variables.data.ellipse_t, I(:,:,2), meanIn);
            I(:,:,3) = ellipseMatrix(variables.data.ellipse_y, variables.data.ellipse_x, variables.data.ellipse_a, variables.data.ellipse_b, variables.data.ellipse_t, I(:,:,3), meanIn);

        end
    end

end