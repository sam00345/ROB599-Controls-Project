clear all
close all
clc
%% Obtain and plot test track
load('TestTrack')
h1=figure;
plot(TestTrack.bl(1,:),TestTrack.bl(2,:),'r')
hold on
plot(TestTrack.br(1,:),TestTrack.br(2,:),'r')
hold on
plot(TestTrack.cline(1,:),TestTrack.cline(2,:),'--b')
hold on
%% Use test track centerline (cline) as desired output specifying where the car needs to be at
ind_section=[1 10];%Index of starting and ending indices of track section to be covered
track_section_x=[TestTrack.bl(1,ind_section); TestTrack.cline(1,ind_section); TestTrack.br(1,ind_section)];
track_section_y=[TestTrack.bl(1,ind_section); TestTrack.cline(1,ind_section); TestTrack.br(1,ind_section)];

%% Decision Variables
%dec vars=[x(k) ;x(k+1) ...;x(k+PredHorizon); u(k); u(k+1);... u(k+PredHorizon-1)];
PredHorizon = 10;
n = 6;%size( A, 1 ); %size of state vector
m = 2;%size( B, 2); % size of input vector
Ndec=n * (PredHorizon+1) + m *PredHorizon ;
%% Inequality constraints
Aineq=zeros(2*m*PredHorizon,Ndec);

%Inputs constraints (remain constant with time)
%Direct less than constraints
Aineq([1:1:m*PredHorizon],n*(PredHorizon+1)+[1:m*PredHorizon])=eye(m*PredHorizon);
bineq([1:2:m*PredHorizon],1)=0.5*ones(m* (PredHorizon)/2,1);
bineq([2:2:m*PredHorizon],1)=6000*ones(m* (PredHorizon)/2,1);

%Greater than constraints formulated as less than consraints
Aineq(m*PredHorizon+[1:1:m*PredHorizon],n*(PredHorizon+1)+[1:m*PredHorizon])=-eye(m*PredHorizon);
bineq(m*PredHorizon+[1:2:m*PredHorizon],1)=0.5*ones(m* (PredHorizon)/2,1);
bineq(m*PredHorizon+[2:2:m*PredHorizon],1)=10000*ones(m* (PredHorizon)/2,1);


% % Cost function
% H = zeros( Ndec,Ndec );
% c = zeros(  Ndec, 1 );
% for k = 1:PredHorizon+1
%     H((k - 1) * n + [1:n], (k - 1) * n + [1:n] ) = Q;
% end
% for i = 1:( PredHorizon - 1)
%     H( n * (PredHorizon+1) + [1:m],n*(PredHorizon+1) + [1:m] ) = R;
% end

%Initial Track boundary constraintsapproximated based on initial velocity
% 
% %Index of closest point on track cline to car x positions
% [mindist,ind]=min(((x(1)-TestTrack.cline(1,:)).^2)+((x(3)-TestTrack.cline(2,:)).^2))
% 
% ux=x(2)*cos(x(5))-x(4)*sin(x(5));
% uy=x(2)*sin(x(5))+x(4)*cos(x(5));
% 
% %Expected x and y locations at next PredHorizon steps
% for kk=1:PredHorizon+1
%   xest(kk)=x(1)+ux*(kk-1)*dt;
%   uest(kk)=x(1)+ux*(kk-1)*dt;
% end
%     
%     
%     if (TestTrack.theta(ind)<pi/2)%If track curves up and right (in 1st quadrant)
%         y_bound_l=interp1(TestTrack.bl(1,:),TestTrack.bl(2,:),dec(1:n:n*(PredHorizon+1)));
%         y_bound_u=interp1(TestTrack.br(1,:),TestTrack.br(2,:),dec(1:n:n*(PredHorizon+1)));
%     else%If track curves up and left (in second quadrant)
%         y_bound_l=interp1(TestTrack.bl(1,:),TestTrack.bl(2,:),dec(1:n:n*(PredHorizon+1)));
%         y_bound_u=interp1(TestTrack.br(1,:),TestTrack.br(2,:),dec(1:n:n*(PredHorizon+1)));
%     end
%     
%     %Direct less than constraints
%     Aineq(2*m*PredHorizon+[1:1:(PredHorizon+1)],[3:n:n*(PredHorizon+1)])=ones(PredHorizon+1,1);
%     bineq(2*m*PredHorizon+[1:1:(PredHorizon+1)],1)=y_bound_l;
%     
%     %Greater than constraints formulated as less than consraints
%     Aineq(2*m*PredHorizon+(PredHorizon+1)+[1:1:(PredHorizon+1)],[3:n:n*(PredHorizon+1)])=-ones(PredHorizon+1,1);
%     bineq(2*m*PredHorizon+(PredHorizon+1)+[1:1:(PredHorizon+1)],1)=-y_bound_u;

%Initial Parameters
x = [287 5 -176 0 2 0]';
u = [0 0]';
dt=0.01;
T=[0:dt:10];
x_all=x;
u_all=[];
dec0=[repmat(x,PredHorizon+1,1);repmat(u,PredHorizon,1)];
UB=zeros(Ndec,1);
UB(1:n*(PredHorizon+1))=100000*ones(n*(PredHorizon+1),1);
UB(n*(PredHorizon+1)+[1:2:m*PredHorizon])=0.5*ones(m*PredHorizon/2,1);
UB(n*(PredHorizon+1)+[2:2:m*PredHorizon])=6000*ones(m*PredHorizon/2,1);

LB=zeros(Ndec,1);
LB(1:n*(PredHorizon+1))=100000*ones(n*(PredHorizon+1),1);
LB(n*(PredHorizon+1)+[1:2:m*PredHorizon])=-0.5*ones(m*PredHorizon/2,1);
LB(n*(PredHorizon+1)+[2:2:m*PredHorizon])=-10000*ones(m*PredHorizon/2,1);

for i = 1:(length(T)-1)
        i
        %Evaluate Al and Bl at current state
        [Al,Bl]=linearized_mats(x',u');
        
        %Discretize (Euler)
        A = eye(size(Al))+dt*Al;
        B = dt*Bl;
        
        %Generate Equality Constraint matrices
        [Aeq, beq] = eq_cons(Al, Bl, x, u, PredHorizon);
                
        nlcons=@(var)track_nlcons(var,TestTrack,Ndec,PredHorizon,n,m);
        fun=@(var)func_cost(var,Ndec,PredHorizon,n,m);
        options=optimoptions(@fmincon,'GradObj','on','MaxIter',1500,'Display','Iter');%'Algorithm','sqp'
        if i>1
        dec0= [reshape(Y(:,2:end),n*(PredHorizon-1),1);Y(:,end);Y(:,end);reshape(u_horiz(:,2:end),m*(PredHorizon-1),1);u_horiz(:,end)];
        end
        
        %Run minimizer
        dec = fmincon(fun,dec0,[],[],Aeq,beq,LB,UB,nlcons,options);
        
        u=dec(n*(PredHorizon+1)+[1:m]);
        u_horiz=reshape(dec(n*(PredHorizon+1)+[1:m*PredHorizon]),m,PredHorizon);
        u_all=[u_all u];
        
        %Forard integrate nonlinear dynamics with input 
        [Y]=nldynamics(u_horiz',x);
        Y=Y';
        x=Y(:,2);
        x_all=[x_all x];
        figure(1)
        hold on
        plot(x(1),x(3),'or');
end
   

save('Results_10s.mat')