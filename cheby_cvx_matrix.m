function [x, R] = cheby_cvx_matrix(A,b)
%The goal is to find the largest Euclidean ball (i.e. its center and
% radius) that lies in a polyhedron described by linear inequalites in this
% fashion: P = {x : a_i'*x <= b_i, i=1,...,m} where x is in R^2
% input:
% A: matrix with 2 rows whose columns are normal vectors of hyperplanes
% b: vecotrs whose entries are shifts of hyperplanes
[m,n]=size(A);
%check the dimension
if m ~= 2
    disp('incorrect dimension')
    return 
end
%compute the norm of each column
A_row=vecnorm(A);
% linear pprogramming in matrix
cvx_begin
    variable r(1)
    variable x_c(2)
    maximize ( r )
    A'*x_c+r*A_row'<= b;
cvx_end
%plot the figure
figure
x = linspace(-2,2);
theta = 0:pi/100:2*pi;
plot( x_c(1) + r*cos(theta), x_c(2) + r*sin(theta), 'r');
hold on
plot(x_c(1),x_c(2),'k+')
%plot each line
for i=1:n
    a=A(:,i);
    c=b(i);
    plot( x, -x*a(1)./a(2) + c./a(2),'b-');
end
xlabel('x_1')
ylabel('x_2')
title('Largest Euclidean ball lying in a 2D polyhedron (matrix)');
axis([-1 1 -1 1])
axis equal
x=x_c;
R=r;
end
