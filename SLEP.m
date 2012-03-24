clear, clc;

% This is an example for running the function glLogisticR
%
%  Problem:
%
%  min  f(x,c) = - weight_i * log (p_i) + rho * sum_j ||x^j||_q
%
%  a_i denotes a training sample,
%      and a_i' corresponds to the i-th row of the data matrix A
%
%  y_i (either 1 or -1) is the response
%     
%  p_i= 1/ (1+ exp(-y_i (x' * a_i + c) ) ) denotes the probability
%
%  weight_i denotes the weight for the i-th sample
%
%  x is grouped into k groups according to opts.ind.
%      The indices of x_j in x is (ind(j)+1):ind(j+1).
%
% For detailed description of the function, please refer to the Manual.
%
% Last modified on August 10, 2009.
%
% For any problem, please contact Jun Liu (j.liu@asu.edu)


% change to the original folder
current_path='/Users/emeryyi/Documents/Research/package/SLEP_package_4.1/';
addpath(genpath([current_path '/SLEP']));
 

m=100;  n=2000;     % The data matrix is of size m x n

ind=[0 100:100:n];   % the indices for the groups
k=length(ind)-1;     % number of groups
q=2;                 % the value of q in the L1/Lq regularization
rho=[0.01:0.01:1];             % the regularization parameter
u=rho(end:-1:1);

randNum=1;           % a random number

% ---------------------- generate random data ----------------------
randn('state',(randNum-1)*3+1);
A=randn(m,n);         % the data matrix

randn('state',(randNum-1)*3+3);
y=[ones(m/2,1);...
    -ones(m/2, 1)];  % the response

%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations

% Normalization
opts.nFlag=0;       % without normalization

% Regularization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)

% Group Property
opts.ind=ind;       % set the group indices
opts.q=q;           % set the value for q
opts.sWeight=[1,1]; % set the weight for positive and negative samples
opts.gWeight=ones(k,1);
                   % set the weight for the group, a cloumn vector

%----------------------- Run the code glLogisticR -----------------------
fprintf('\n mFlag=0, lFlag=0 \n');
opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search
opts.fName='glLogisticR';
tic;
%[x1, c1, funVal1, ValueL1]= glLogisticR(A, y, rho, opts);
[X, C]=pathSolutionLogistic(A, y, u, opts);
toc;
