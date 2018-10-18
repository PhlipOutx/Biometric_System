rho=0.90; % Use correlation coefficient rho=0.9 as default
STFTalpha=1/winSize;  % alpha in STFT approaches (for Gaussian derivative alpha=1) 
r=(winSize-1)/2; % Get radius from window size
x=-r:r; % Form spatial coordinates in window
u=1:r; % Form coordinates of positive half of the Frequency domain (Needed for Gaussian derivative)
convmode='valid';

w0=(x*0+1);
w1=exp(complex(0,-2*pi*x*STFTalpha)); % 0 real component, other part is imaginary % exp is exponential for each element
w2=conj(w1); % complex conjugate
	
%% Run filters to compute the frequency response in the four points. Store real and imaginary parts separately
% Run first filter
filterResp=conv2(conv2(img,w0.',convmode),w1,convmode);
% Initilize frequency domain matrix for four frequency coordinates (real and imaginary parts for each frequency).
freqResp=zeros(size(filterResp,1),size(filterResp,2),8); 
% Store filter outputs
freqResp(:,:,1)=real(filterResp);
freqResp(:,:,2)=imag(filterResp);
% Repeat the procedure for other frequencies
filterResp=conv2(conv2(img,w1.',convmode),w0,convmode);
freqResp(:,:,3)=real(filterResp);
freqResp(:,:,4)=imag(filterResp);
filterResp=conv2(conv2(img,w1.',convmode),w1,convmode);
freqResp(:,:,5)=real(filterResp);
freqResp(:,:,6)=imag(filterResp);
filterResp=conv2(conv2(img,w1.',convmode),w2,convmode);
freqResp(:,:,7)=real(filterResp);
freqResp(:,:,8)=imag(filterResp);

% Read the size of frequency matrix
[freqRow,freqCol,freqNum]=size(freqResp);




%% If decorrelation is used, compute covariance matrix and corresponding whitening transform
if decorr == 1
    % Compute covariance matrix (covariance between pixel positions x_i and x_j is rho^||x_i-x_j||)
    [xp,yp]=meshgrid(1:winSize,1:winSize);
    pp=[xp(:) yp(:)];
    dd=dist(pp,pp');
    C=rho.^dd;
    
    % Form 2-D filters q1, q2, q3, q4 and corresponding 2-D matrix operator M (separating real and imaginary parts)
    q1=w0.'*w1;
    q2=w1.'*w0;
    q3=w1.'*w1;
    q4=w1.'*w2;
    u1=real(q1); u2=imag(q1);
    u3=real(q2); u4=imag(q2);
    u5=real(q3); u6=imag(q3);
    u7=real(q4); u8=imag(q4);
    M=[u1(:)';u2(:)';u3(:)';u4(:)';u5(:)';u6(:)';u7(:)';u8(:)'];
    
    % Compute whitening transformation matrix V
    D=M*C*M';
    A=diag([1.000007 1.000006 1.000005 1.000004 1.000003 1.000002 1.000001 1]); % Use "random" (almost unit) diagonal matrix to avoid multiple eigenvalues.  
    [U,S,V]=svd(A*D*A);
   
    % Reshape frequency response
    freqResp=reshape(freqResp,[freqRow*freqCol,freqNum]);

    % Perform whitening transform
    freqResp=(V.'*freqResp.').';
    
    % Undo reshape
    freqResp=reshape(freqResp,[freqRow,freqCol,freqNum]);
end



%% Perform quantization and compute LPQ codewords
LPQdesc=zeros(freqRow,freqCol); % Initialize LPQ code word image (size depends whether valid or same area is used)
for i=1:freqNum
    LPQdesc=LPQdesc+(double(freqResp(:,:,i))>0)*(2^(i-1));
end

%% Histogram if needed
if strcmp(mode,'nh') || strcmp(mode,'h')
    LPQdesc=hist(LPQdesc(:),0:255);
end
%% Normalize histogram if needed
if strcmp(mode,'nh')
    LPQdesc=LPQdesc/sum(LPQdesc);
end