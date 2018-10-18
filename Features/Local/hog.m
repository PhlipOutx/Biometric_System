function hogd = hog(in_img, nbins)
%%%% function to compute gradient orientaion histogram of an image
%%%% (grayscale or color). It returns a pyramidal descriptor, the details
%%%% are specfied in the function.

	%nbins = 12; %% number of histogram bins 
	L = 0; %% pyramid levels
	rscale_dim = 25; %% rescaling dimension for input images
	mdfilt_ker = [5,5]; %% median filter kernel
	buf = 1; %% boundary to exclude (in pixel width)
	ulight = 0; %% if 1 - perform gamma correction in case of uncontrolled lighting
	gamma = 1.5; %% gamma correction value

	%load eye_mask; %% if masking of the eye is necessary, uncomment this

	[nrows1, ncols1, ncha] = size(in_img);
	if ncha > 1
		J = double(rgb2gray(in_img));
	else
		J = double(in_img);
	end
	%Jr = RescaleImage(J, rscale_dim, rscale_dim); %% rescale image
	Jr = J;
	%mask = RemoveSpecularReflections(Jr); %% preprocess
	%mask = 1 - mask;
	mask = 1;

	%mask = mask & eye_mask; %% uncomment if the eye has to be masked

	I = medfilt2(Jr, mdfilt_ker); %% smooth by median filtering (uncomment if not helpful)
	% if ulight ~= 1
	%     I = imadjust(uint8(Jr),stretchlim(uint8(Jr)),[]);
	% else
	%     Jr = (Jr - min(min(Jr)))/(max(max(Jr)) - min(min(Jr)));
	%     I = Jr.^(1/gamma);
	%     I = uint8(I.*255);
	% end
	[nrows, ncols, ncha] = size(I);
	%[gmag, gang] = ApplySteerableFilters(double(I), nrows, ncols, 2.0);
	ker = fspecial('prewitt');
	gx = conv2(double(I), ker, 'same');
	gy = conv2(double(I), ker', 'same');
	gmag = sqrt((gx.*gx) + (gy.*gy));
	gang = atan2(gy, gx);

	roi = [1+buf, nrows - buf, 1+buf, ncols - buf];

	gmag(1:buf,:) = min(min(gmag));
	gmag(:, 1:buf) = min(min(gmag));
	gmag(nrows-(1-buf):nrows, :) = min(min(gmag));
	gmag(:, ncols-(1-buf):ncols) = min(min(gmag));
	gmag = (gmag - min(min(gmag)))/(max(max(gmag)) - min(min(gmag)));

	max_ang = max(max(gang));
	[rr, cc] = find(gang < 0);
	for i = 1:length(rr)
		gang(rr(i,1), cc(i,1)) = max_ang + (pi - abs(gang(rr(i,1), cc(i,1))));
	end

	bin_id_map = GetWeightedOrientation(gmag, gang, mask, nbins);

	hogd = ComputePhogDescriptor(bin_id_map(roi(1,1):roi(1,2), roi(1,3):roi(1,4)), gmag(roi(1,1):roi(1,2), roi(1,3):roi(1,4)), L, nbins);
	
end

function bin_id_map = GetWeightedOrientation(gmag, gang, mask, nbins)

	binw = 2*pi/nbins;
	bin_id_map = gang./binw;
	bin_id_map = floor(bin_id_map) + 1;

end

function p = ComputePhogDescriptor(bh,bv,L,bin)
% anna_PHOGDESCRIPTOR Computes Pyramid Histogram of Oriented Gradient over a ROI.
%               
%IN:
%	bh - matrix of bin histogram values
%	bv - matrix of gradient values 
%   L - number of pyramid levels
%   bin - number of bins
%
%OUT:
%	p - pyramid histogram of oriented gradients (phog descriptor)

	p = [];
	%level 0
	for b=1:bin
		ind = bh==b;
		p = [p;sum(bv(ind))];
	end
			
	cella = 1;
	for l=1:L
		x = fix(size(bh,2)/(2^l));
		y = fix(size(bh,1)/(2^l));
		xx=0;
		yy=0;
		while xx+x<=size(bh,2)
			while yy +y <=size(bh,1) 
				bh_cella = [];
				bv_cella = [];
				
				bh_cella = bh(yy+1:yy+y,xx+1:xx+x);
				bv_cella = bv(yy+1:yy+y,xx+1:xx+x);
				
				for b=1:bin
					ind = bh_cella==b;
					p = [p;sum(bv_cella(ind))];
				end 
				yy = yy+y;
			end        
			cella = cella+1;
			yy = 0;
			xx = xx+x;
		end
	end
	if sum(p)~=0
		p = p/sum(p);
	end
end

