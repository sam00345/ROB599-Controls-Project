function main
% optimal control control using the collocation approach
close all; clear all;
global N tau0 tauf t0 tf nx nu
t0   = 0;    % true initial time
tau0 = 0;    % scaled initial time
tauf = 1;    % scaled final time
N    = 20;   % number of collocation nodes
h    = (tauf - tau0)/(N - 1);

nx = 3; % number of states
nu = 1; % number of controls

par0(1:nx*N,1) = 0;
par0 = [par0;zeros(nu*N,1);tauf];

options = optimset('TolX',1e-8,'Algorithm','sqp','Display','iter',...
                            'MaxFunEvals',10000,'MaxIter', 1000);

[paropt,fval,exitflag,output,lambda] = fmincon(@col_cost,par0,...
                             [],[],[],[],[],[],@col_con,options);

[cost, XUTraj] = col_cost(paropt);

% plot
tf = paropt(end);
figure(1);plot((tau0:h:tauf)*(tf-t0),XUTraj(:,4),'bo-', 'linewidth',3);
xlabel('Time (sec)','fontsize',16);ylabel('u','fontsize',16);
set(gca,'fontsize',16)
figure(2); quiver(XUTraj(:,1),XUTraj(:,2),cos(XUTraj(:,3)),sin(XUTraj(:,3)),'g');
hold on;
plot(XUTraj(:,1),XUTraj(:,2),'bo-','linewidth',3);
xlabel('x_1','fontsize',16);ylabel('x_2','fontsize',16);
set(gca,'fontsize',16)

return;

function [cost,XUTraj] = col_cost(par)
% cost function - used with the collocation method

global N nx nu t0 tauf tau0

X    = par(1:nx*N);
U    = par(nx*N+1:end-1);
tf   = par(end);
h    = (tauf - tau0)/(N - 1);
cost = 0;
for (i = 1:1:N),
    ui    = U(nu*(i-1)+1:nu*i);
    xi    = X(nx*(i-1)+1:nx*i);
    cost  = cost + 1/2*norm(ui,2)^2*(tf-t0)*h; % 1/2*sqrt(ui^2+1e-6)*(tf-t0)*h; 
    XUTraj(i,:) = [xi(:)',ui];
end;
return;

function [Ci,Ce] = col_con(par)
% constraint function - used with the collocation method

global N tauf tau0 nx nu t0

Ci   = [];
Ce   = [];

h    = (tauf - tau0)/(N - 1);
X    = par(1:nx*N);
U    = par(nx*N+1:end-1);
tf   = par(end);

for (i=1:1:N-1),
    
    xi    = par(nx*(i-1)+1:nx*i);
    xip1  = par(nx*(i)+1:nx*(i+1));
    
    ui    = U(nu*(i-1)+1:nu*i);
    uip1  = U(nu*i+1:nu*(i+1));
    
    fi     = odemodel_col(xi,ui)*(tf - t0);
    fip1   = odemodel_col(xip1,uip1)*(tf - t0);
    
    fistar = odemodel_col(1/2*(xi + xip1)+h/8*(fi - fip1), 1/2*(ui + uip1))*(tf - t0);
    
    Ce(nx*(i-1)+1:nx*i) = -3/2*(xi - xip1) - h/4*(fi + fip1) - h*fistar;
    
end;

x0 = [0; 0; pi/2]; % specified initial state
xf = [5; 0; pi/2]; % desired final state
Ce = Ce(:);
Ce = [-Ce;...
    -x0 + par(1:3);...
    par(3*(N-1)+1) - xf(1);...
    par(3*(N-1)+2) - xf(2);...
    par(3*(N-1)+3) - xf(3)];
Ci = [-(tf - t0)];
return;

function xdot = odemodel_col(x,u)
% ODE model (in actual/non-scaled time) used with collocation method

x1        = x(1);
x2        = x(2);
x3        = x(3);

x1dot     = cos(x3);
x2dot     = sin(x3);
x3dot     = u;
xdot      = [x1dot,x2dot,x3dot]';
return;