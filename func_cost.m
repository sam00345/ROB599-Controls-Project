function [cost,J] = func_cost(x,Ndec,PredHorizon,n,m)
% Maximize (Minimize the negative of) distance travelled in given time
cost=0;
J=zeros(Ndec,1);
for kk=2:(PredHorizon+1)
cost=cost-(x((kk)*n+1)-x((kk-1)*n+1))^2-(x((kk)*n+3)-x((kk-1)*n+3))^2;
dJ=zeros(Ndec,1);
dJ((kk)*n+1,1)=-2*(x((kk)*n+1)-x((kk-1)*n+1));
dJ((kk-1)*n+1,1)=2*(x((kk)*n+1)-x((kk-1)*n+1));
dJ((kk)*n+3,1)=-2*(x((kk)*n+3)-x((kk-1)*n+3));
dJ((kk-1)*n+3,1)=2*(x((kk)*n+3)-x((kk-1)*n+3));
J=J+dJ;
end

