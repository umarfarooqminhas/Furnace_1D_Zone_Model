
% Thomas Algorithm for tridiagonal system of equations of the form [a][x] = [b] % By SR
% n - length of vector b/rows/columns of a square matrix a
%
% d - diagnoal vector of the matrix A
% u- upper diagnoal vector of the matrix A
%1 - lower diagnoal vector of the matrix A % b - rhs vector
function x = Thomas_Algorithm_Func1(d,l,u,b)
n = length(b);
u(n) = 0;
d = d';
l=l';
u = u';
%new = [l d u];
for i = 2:n
    d(i) = d(i) - (l(i)/d(i-1))*u(i-1);
    b(i) = b(i) - (l(i)/d(i-1))*b(i-1);
end

x(n) = b(n)/d(n);
for k = 1:n-1
i=n-k;
x(i) = (b(i) - u(i)*x(i+1))/d(i);
end