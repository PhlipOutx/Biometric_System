% The function computes features for nonlinear subspace projection techniques, e.g., KPCA, KFA, etc.
% 
% PROTOTYPE
% feat = nonlinear_subspace_projection_PhD(X, model);
% 
% USAGE EXAMPLE(S)
% 
%     Example 1:
%       %we generate random training data for the example
%       X=double(imread('sample_face.bmp'));
%       training_data = randn(length(X(:)),10);
%       model = perform_kpca_PhD(training_data, 'poly', [0 2]);
%       feat = nonlinear_subspace_projection_PhD(X(:), model);
%
%
% GENERAL DESCRIPTION
% The function computes features (or subspace projection coefficients)
% for kernel (non-linear) subspace projection techniques, such as
% KPCA, KFA, and others, given a model containing the parameters of the 
% subspace method. Note that this function does not compute the transforms 
% (i.e., the subspaces, transformation matrices, etc.), it only uses a 
% precomputed model to calculate the subspace projection.
% 
% Here, the input (test) data in X has to be arranged into a column data matrix.
% This means that if the test-data matrix X is of size BxA, then there are
% A test samples (images) each having B variables (pixels). Note, however,
% that the dimensionality of the input data has to match that of the
% training data used to produce the kernel subspace.
% 
% The input data can be a single image (in vector) form, in which case the
% the output "feat" contains the feature vector corresponding to the input 
% test sample X, or it can be a data matrix, in which case the output "feat" 
% is also a matrix, containing in its columns the feature vectors 
% corresponding to the test samples in X. 
%
% To see an example of usage of the function look at the demo folder.
%
% 
% REFERENCES
% The function implements a simple projection of the test data into a
% low-dimensioanl subspace, which is documented in a lot of references, see
% for example:
% 
% Štruc V., Pavešic, N.:Gabor-Based Kernel Partial-Least-Squares 
% Discrimination Features for Face Recognition, Informatica (Vilnius), vol.
% 20, no. 1, pp. 115-138, 2009.
%
%
%
% INPUTS:
% X                     - test-data matrix of size BxA, where each of
%                         the columns contains A samples, each having a
%                         dimensionality of B
% model                 - the subspace model generated using the
%                         appropriate function (see the RELATED FUNCTIONS 
%                         part of this help)
%
% OUTPUTS:
% feat                  - a vector or matrix (depending on the input)
%                         containing the subspace features of the input 
%                         test samples 
%
% NOTES / COMMENTS
% The dimensionality of the feature vectors produced here is defined by the
% dimensionality of the transformation matrix in "model.W". I have not
% included an option so far to reduce the dimensionality of the feature 
% vectors at this point. I will do that in the classification stage. 
%
% The function was tested with Matlab ver. 7.9.0.529 (R2009b) and Matlab 
% ver. 7.11.0.584 (R2010b).
% 
% 
% RELATED FUNCTIONS (SEE ALSO) 
% perform_kpca_PhD 
% perform_kfa_PhD
% linear_subspace_projection_PhD
% 
% 
% ABOUT
% Created:        15.2.2010
% Last Update:    24.11.2011
% Revision:       1.0
% 
%
% WHEN PUBLISHING A PAPER AS A RESULT OF RESEARCH CONDUCTED BY USING THIS CODE
% OR ANY PART OF IT, MAKE A REFERENCE TO THE FOLLOWING PUBLICATIONS:
% 
% Štruc V., Pavešic, N.: The Complete Gabor-Fisher Classifier for Robust 
% Face Recognition, EURASIP Advances in Signal Processing, vol. 2010, 26
% pages, doi:10.1155/2010/847680, 2010.
%
% Štruc V., Pavešic, N.:Gabor-Based Kernel Partial-Least-Squares 
% Discrimination Features for Face Recognition, Informatica (Vilnius), vol.
% 20, no. 1, pp. 115-138, 2009.
% 
% 
% The BibTex entries for the papers are here
% 
% @Article{ACKNOWL1,
%     author = "Vitomir \v{S}truc and Nikola Pave\v{s}i\'{c}",
%     title  = "The Complete Gabor-Fisher Classifier for Robust Face Recognition",
%     journal = "EURASIP Advances in Signal Processing",
%     volume = "2010",
%     pages = "26",
%     year = "2010",
% }
% 
% @Article{ACKNOWL2,
%     author = "Vitomir \v{S}truc and Nikola Pave\v{s}i\'{c}",
%     title  = "Gabor-Based Kernel Partial-Least-Squares Discrimination Features for Face Recognition",
%     journal = "Informatica (Vilnius)",
%     volume = "20",
%     number = "1",
%     pages = "115–138",
%     year = "2009",
% }
% 
% Official website:
% If you have down-loaded the toolbox from any other location than the
% official website, plese check the following link to make sure that you
% have the most recent version:
% 
% http://luks.fe.uni-lj.si/sl/osebje/vitomir/face_tools/PhDface/index.html
%
% 
% OTHER TOOLBOXES 
% If you are interested in face recognition you are invited to have a look
% at the INface toolbox as well. It contains implementations of several
% state-of-the-art photometric normalization techniques that can further 
% improve the face recognition performance, especcially in difficult 
% illumination conditions. The toolbox is available from:
% 
% http://luks.fe.uni-lj.si/sl/osebje/vitomir/face_tools/INFace/index.html
% 
%
% Copyright (c) 2011 Vitomir Štruc
% Faculty of Electrical Engineering,
% University of Ljubljana, Slovenia
% http://luks.fe.uni-lj.si/en/staff/vitomir/index.html
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files, to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
% 
% November 2011

function features = nonlinear_subspace_projection(X, model);
	fprintf('Apply linear subspace projection\n');
	%%% Initilize variables %%%
	features=[];
	
	% check the size of the test data
	[a,b]=size(X);
	[c,d]=size(model.X);
	assert(c == a, 'The dimensionality of the test data does not match that of the data used for training');

	%%% perform subspace projection %%%
	[c,d]=size(model.J);

	% compute the test kernel matrix
	fprintf('1. Compute test data kernel matrix\n');
	t3 = tic;
	K = compute_kernel_matrix(X,model.X,model.typ,model.args);
	toc(t3)

	% center the test kernel matrix
	fprintf('2. Centering data\n');
	t3 = tic;
	Jt=ones(b,1)/c;
	J=ones(c,1);
	K=K-((Jt*J')*model.K)-(K*(1/c)*(J*J'))+(Jt*J'*model.K*(1/c)*(J*J'));
	K=K';
	toc(t3)

	% feature computation
	fprintf('3. Applying subspace\n');
	t3 = tic;
	features = model.W'*K;
	toc(t3)

end

