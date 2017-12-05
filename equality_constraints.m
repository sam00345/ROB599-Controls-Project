%Some test parameters
x = [287 5 -176 0 2 0]';
u = [0 0]';
horizon = 10;
k = 1;
dt = 0.01;

%Evaluate Al and Bl at current state
[Al,Bl]=linearized_mats(x(:,k)',u(:,k)');

%Discretize (Euler)
A = eye(size(Al))+dt*Al;
B = dt*Bl;

%Generate Equality Constraint matrices
[Aeq, beq] = eq_cons(Al, Bl, x(:,k), u(:,k), horizon);

function [Aeq, beq] = eq_cons(A, B, x_k, u_k, horizon)

n = length(x_k);                    %Number of States
m = length(u_k);                    %Number of Inputs

xsize = (horizon+1)*n;              %Number of State Decision Variables
zsize = (horizon+1)*n+horizon*m;    %Total Number of Decision Variables

Aeq = zeros(xsize, zsize);          %Allocate Aeq
Aeq(1:n, 1:n) = eye(n);             %Initial Condition LHS

beq = zeros(xsize, 1);              %Allocate beq
beq(1:n) = x_k;                     %Initial Condition RHS

j = xsize+1;

for i = n+1:n:xsize
    
    Aeq(i:i+n-1, i:i+n-1) = -eye(n);    %x(k+1) term
    Aeq(i:i+n-1, i-n:i-1) = A;          %A*x(k) term
    Aeq(i:i+n-1, j:j+m-1) = B;          %B*u(k) term
    
    j = j+m;
    
end

end