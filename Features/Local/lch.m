function result=lch(varargin)%pass in image, get back normalized histogram
% Check number of input arguments.
error(nargchk(1,3,nargin));

image=varargin{1};
d_image=double(image);

if nargin==1
    dim = [1 2];
    bounds = [0 1; 0 1];
end

if (nargin == 3)
    dim = varargin{2};
    bounds = varargin{3};
end

result=colorHist2d(d_image, dim, bounds);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hist]=colorHist2d(img, dim, bounds)
%creates a 2d histogram of an image split into 4x4 bins
%returns flat histogram vector of size 16
%img is the true color image, should be a double with range 0-1
%dim is a 1x2 vector of the dimensions of the image used to make the histogram
%bounds is a 2x2 [lbound ubound; lbound2 ubound2]
% usage & matlab bounds
%   colorHist(Ihsv,[1 2],[0 1; 0 1])---v [0 1]
%   colorHist(Iycbcr,[2 3],[16/255 240/255; 16/255 240/255])---y [16/255 235/255]
%   colorHist(Irgb,[1 2],[0 255;0 255]) or colorHist(Irgb,[1 2],[0 1; 0 1])
%   colorHist(Iyiq,[2 3],[-.5959 .5959; -.5229 .5229]) --- y [0 1]
%   colorHist(Iluv,[2 3],[-100 100; -100 100]);--- l [0 100]
numBins=4;
[r c]=size(img(:,:,1));
diff1=bounds(1,2)-bounds(1,1);
diff2=bounds(2,2)-bounds(2,1);
step1=diff1/numBins;
step2=diff2/numBins;

bins(1,:)=bounds(1,1):step1:bounds(1,2);
bins(2,:)=bounds(2,1):step2:bounds(2,2);

hist=zeros(1,numBins*numBins);

for i=1:r
    for j=1:c
        ndx=zeros(1,2);
        val=zeros(1,2);
        for d=1:2
            val(d)=img(i,j,dim(d));

            if ((val(d)>=bins(d,1) && val(d)<bins(d,2)) )
                ndx(d)=1;
            elseif(val(d)>=bins(d,2) && val(d)<bins(d,3))
                ndx(d)=2;
            elseif(val(d)>=bins(d,3) && val(d)<bins(d,4))
                ndx(d)=3;
            elseif(val(d)>=bins(d,4) && val(d)<=bins(d,5))
                ndx(d)=4;
            end
        end
        if(ndx(1)>0 && ndx(2)>0)
            off1=numBins*(ndx(1)-1);
            off2=ndx(2);
            off=off1+off2;
            hist(off)=hist(off)+1;
        end
    end
end

%hist=hist/(r*c);


