function WLDdesc = wld(img,num_orientation,num_excitation)
% Funtion WLDdesc=wld(img,num_orientation,num_excitation) computes the Weber Local Descriptor (WLD)
% for the input image img. Descriptors are calculated using only valid pixels.
%
% Inputs:
% img = N*N uint8 or double, format gray scale image to be analyzed.
% num_orientation = 1*1 double, number of orientation levels.
% num_excitation = 1*1 double, number of excitation levels.
%
% Output:
% WLDdesc = 1*num_orientation*num_excitation double or uint8, WLD descriptors histogram


%% Default parameters
alpha = 6;
beta = 0;
WLDdesc = zeros(num_excitation,num_orientation);

for i=2:size(img,1)-1
	for j=2:size(img,2)-1
		
		% Differential Excitation
		sigma = atan2(double((img(i-1,j-1)+img(i-1,j)+img(i-1,j+1)+img(i,j+1)+img(i+1,j+1)+img(i+1,j)+img(i+1,j-1)+img(i,j-1)-(img(i,j)*8))*alpha),double(img(i,j)+beta));
		sigma = radtodeg(sigma);
		%if ( sigma > 180.0 )
		%	sigma = sigma - 360.0;
		%end
		c = ceil( (sigma+0.001) / 360.0 * num_excitation );
		
		% Orientation
		theta = atan2(double(img(i+1,j) - img(i-1,j)), double(img(i,j-1) - img(i,j+1)));
		theta = radtodeg(theta);
		t = ceil( (theta+0.001) / 360.0 * num_orientation );
		
		WLDdesc(c,t) = WLDdesc(c,t)+1;
	end
end
WLDdesc = WLDdesc(:);
