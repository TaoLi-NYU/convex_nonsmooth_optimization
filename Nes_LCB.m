function [f,g] = Nes_LCB(x,m,M,T,e1)
%Nesterove's worst case example
%input: m,M the lower and upper bound of eigenvalues; T: the tri-diagnoal matrix in gradient and heesian   
%compute the function evaluations
main=norm(x(1:end-1)-x(2:end))^2;
main=x(1)^2+main-2*x(1);
f=(M-m)/8*main+m/2*norm(x)^2;
%compute the gradient
g=(M-m)/4*T*x+m*x-(M-m)/4*e1;
end

