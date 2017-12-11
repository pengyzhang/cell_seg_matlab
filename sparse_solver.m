function [x_l1, yssc] = sparse_solver(y, dict, methodIndex)
addpath('.\L1Solvers');
addpath('.\l1_ls_matlab');
alpha=1;
beta=0.17;
lambda=0.01;
[M,N] = size(dict);
D = [dict eye(M)];
if (length(who('methodIndex'))<=0)
    methodIndex=1;
end
if(methodIndex==1)
    [lambda_max] = find_lambdamax_l1_ls(D',y);                   
    lamb         = lambda_max*lambda;
    lamb         = 0.01*lamb;
    [xEst status] = l1_ls(D, y, lamb,1e-3,true);
elseif(methodIndex==2)
    tic;
    [xEst, nIter, curTimeEst, curErrEst] = SolveHomotopy(D, y,N,alpha,beta,'maxiteration', 200) ;
    tEst = toc ; 
end
x_l1=xEst(1:N);
yssc=dict*x_l1;