figure(1)
%#itertions 
iter1=1:length(f_all);
iter2=1:length(f_alln);
%optimal values
p1=f_all(end);
p2=f_alln(end);
%logplot
semilogy(iter1,f_all-p1,'ro');
hold on
semilogy(iter2,f_alln-p2,'go');
legend('gradient method', 'Newton')
title('logplot of fucntion values using different methods')
xlabel('#iterations')
ylabel('log of function values')
%%
figure(2)
semilogy(iter1,gnorm_all,'rx');
hold on
semilogy(iter2,gnorm_alln,'gx');
legend('gradient method', 'Newton')
title('logplot of norms of gradients using different methods')
xlabel('#iterations')
ylabel('log of norm of gradients')