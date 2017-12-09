function [cost,J] = func_cost2(x,Ndec,PredHorizon,n,m)
% Maximize (Minimize the negative of) distance travelled in given time
cost=0;
J=zeros(Ndec,1);
cost=(1472-x(n*PredHorizon+1))^2+(817.8-x(n*PredHorizon+3))^2;
% dJ=zeros(Ndec,1);
J(n*PredHorizon+1,1)=-2*(1472-x(n*PredHorizon+1));
J(n*PredHorizon+3,1)=-2*(817.8-x(n*PredHorizon+3));
