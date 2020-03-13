function [f_all, gnorm_all, h_all,x_all]=newtmeth(fun, x0,tol, maxit)
alpha = 0.25;
beta = 0.5;
count = 0;
e=1;
%pre-allocation
f_all=zeros(maxit,1);
gnorm_all=zeros(maxit,1);
h_all=zeros(length(x0),length(x0),maxit);
x_all=zeros(length(x0),maxit);
% error above the tol 
while count<maxit && e>tol
    count=count+1;
    t=1;
    [f,g,h]=fun(x0);
    f_all(count)=f;
    e=norm(g);
    gnorm_all(count)=e;
    h_all(:,:,count)=h;
    x_all(:,count)=x0;
    %newton method
    delta=-h\g;
    %backline tracking
    while true 
        x1=x0+t*delta;
        [f1,~,~]=fun(x1);
        %armijo condition
        f_arm=f-alpha*t*g'*delta;
        if f1>f_arm
            t=t*beta;
        else
            break;
        end
    end
    x0=x1;
end
%we only keep the valid iterations
f_all=f_all(1:count);
gnorm_all=gnorm_all(1:count);
h_all=h_all(:,:,1:count);
x_all=x_all(:,1:count);
p=f_all(end);
fprintf('the minimal value is %f',p)
iter=1:count;
figure(1)
semilogy(iter,f_all-p,'ro');
title('logplot of fucntion values')
xlabel('#iterations')
ylabel('log of function values')
figure(2)
semilogy(iter,gnorm_all,'bo');
title('logplot of norms of gradients')
xlabel('#iterations')
ylabel('log of norms of gradients')

end