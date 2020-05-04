function [Xcomplete,norm_diff,obj] = admm_complete(X,num_entries,rho, tol, maxit)
% We convert the SDP into a LP problem and use ADMM to solve it. 
%   inputs: X: original matrix; num_entries: number of known entries, rho:
%   penalty parameter; tol: tolerance; maxit: maximum iteration numbers
% dimensions
M=size(X,1);
N=size(X,2);
% construct the vector q in the objective function
%eye_block=[eye(M),zeros(M,N);zeros(N,M),eye(N)];
q=reshape(eye(M+N),[(M+N)^2,1]);
% randomly generate row and col index
rng(1)
row_index=randi(M,num_entries,1);
col_index=randi(N,num_entries,1);
% construct A and b in the programming problem
b=zeros(num_entries,1);
A=zeros(num_entries,(M+N)^2);
for k=1:num_entries
    i=row_index(k);
    j=col_index(k);
    % b(k) is the know entry
    b(k)=X(i,j);
    e_ij=zeros(M+N,M+N);
    e_ij(i,M+j)=0.5;
    e_ij(M+j,i)=0.5;
    A(k,:)=reshape(e_ij,[(M+N)^2,1]);
end
%
A=sparse(A);
q=sparse(q);
I=speye((M+N)^2);
%ADMM
tol_abs=tol(1);
tol_rel=tol(2);
%initialization
count=0;
z0=zeros((M+N)^2,1);
u0=zeros((M+N)^2,1);

B=[rho*I, A';A, zeros(num_entries)];
tic
while count<maxit
% x-update
%solving KKT system for x
    rhs=[rho*(z0-u0)-q; b];
    x1_temp=B\rhs;
    x1=x1_temp(1:(M+N)^2);
% z-update projection
    z1_reshape=reshape(x1+u0,[M+N,M+N]);
    [V,D]=eig(z1_reshape);
    %projection
    D(D<0)=0;
    z1_temp=V*D*V';
    %z1_temp(z1_temp<1e-14)=0;
    z1=reshape(z1_temp,[(M+N)^2,1]);
% u-update 
    u1=u0+x1-z1;
    count=count+1;
% stopping criterion
    %primal residual
    pr=norm(x1-z1);
    %dual residual
    dr=rho*norm(z1-z0);
    e_pri=sqrt(N)*tol_abs+tol_rel*max(norm(x1),norm(z1));
    e_dual=sqrt(M)*tol_abs/rho+tol_rel*norm(u1);
    if pr<e_pri && dr<e_dual
        disp('success')
        break
    end
    %assign new values 
    z0=z1;
    u0=u1;
end
toc
% extract the completed matrix 
xtilde=reshape(x1,[M+N,M+N]);
Xcomplete=xtilde(1:M,M+1:end);
norm_diff=norm(X-Xcomplete,'fro');
obj=q'*x1;
end


    
    
    