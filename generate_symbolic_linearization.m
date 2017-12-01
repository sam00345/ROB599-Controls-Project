clear all 
close all
clc
%%
W=13720;
Nw=2;
f=0.01;
Iz=2667;
a=1.35;
b=1.45;
By=0.27;
Cy=1.2;
Dy=2921;
Ey=-1.6;
Shy=0;
Svy=0;
m=1400;

syms phi_yf phi_yr a_f a_r u v r x y psvar F_x delta_f F_yr

a_f= 180/pi*(delta_f-atan2(v+a*r,u));
a_r= 180/pi*(-atan2((v-b*r),u));

phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));


F_yf= Dy*sin(Cy*atan(By*phi_yf))+Svy;
F_yr= Dy*sin(Cy*atan(By*phi_yr))+Svy; 

f=[u*cos(psvar)-v*sin(psvar);...
          (-f*W+Nw*F_x-F_yf*sin(delta_f))/m+v*r;...
          u*sin(psvar)+v*cos(psvar);...
          (F_yf*cos(delta_f)+F_yr)/m-u*r;...
          r;...
          (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
      
sym('A_lin',[6,6]);
sym('B_lin',[6,2]);

A_lin(1:6,1)=diff(f,x);
A_lin(1:6,2)=diff(f,u)
A_lin(1:6,3)=diff(f,y);
A_lin(1:6,4)=diff(f,v);
A_lin(1:6,5)=diff(f,psvar);
A_lin(1:6,6)=diff(f,r);

B_lin(1:6,1)=diff(f,delta_f);
B_lin(1:6,2)=diff(f,F_x);

save('Linearized_expressions.mat','A_lin','B_lin');


% eval(subs(A_lin,[x u y v psvar r delta_f F_x],[x0 [0 0]]))
% eval(subs(B_lin,[x u y v psvar r delta_f F_x],[x0 [0 0]]))