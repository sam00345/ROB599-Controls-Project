function [c,ceq] = track_nlcons(x,TestTrack,Ndec,PredHorizon,n,m)

orig=[200,-200];%Origin for getting vectors
%Vecor to CM of car must lie inside track vector
c=zeros(2*(PredHorizon+1),1);% Compute nonlinear inequalities at x.
for kk=1:(PredHorizon+1)
    %Find closest point on track
    [mindist,ind]=min(((x((kk-1)*n+1)-TestTrack.cline(1,:)).^2)+((x((kk-1)*n+3)-TestTrack.cline(2,:)).^2));
    c(2*kk-1,1) = atan2((x((kk-1)*n+3)-orig(2)),(x((kk-1)*n+1)-orig(1)))-atan2((TestTrack.bl(2,ind)-orig(2)),(TestTrack.bl(1,ind)-orig(1))); 
    c(2*kk,1) = -atan2((x((kk-1)*n+3)-orig(2)),(x((kk-1)*n+1)-orig(1)))+atan2((TestTrack.br(2,ind)-orig(2)),(TestTrack.br(1,ind)-orig(1))); % Compute nonlinear inequalities at x.
end

ceq = [];   % Compute nonlinear equalities at x.