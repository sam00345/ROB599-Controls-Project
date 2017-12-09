function [c,ceq] = track_nlcons(x,TestTrack,Ndec,PredHorizon,n,m)

orig=[200,-200];%Origin for getting vectors
%Vecor to CM of car must lie inside track vector
c=zeros((PredHorizon+1),1);% Compute nonlinear inequalities at x.
for kk=1:(PredHorizon+1)
    %Find closest point on track
    [mindist,ind]=sort(((x((kk-1)*n+1)-TestTrack.cline(1,:)).^2)+((x((kk-1)*n+3)-TestTrack.cline(2,:)).^2),'ascend');
    %Check if 
    [in] = inpolygon(x((kk-1)*n+1),x((kk-1)*n+3),[TestTrack.bl(1,ind(1)) TestTrack.bl(1,ind(2)) TestTrack.br(1,ind(2)) TestTrack.br(1,ind(1))],[TestTrack.bl(2,ind(1)) TestTrack.bl(2,ind(2)) TestTrack.br(2,ind(2)) TestTrack.br(2,ind(1))]);
    if in==1
    c(kk,1) = -1;%atan2((x((kk-1)*n+3)-orig(2)),(x((kk-1)*n+1)-orig(1)))-atan2((TestTrack.bl(2,ind)-orig(2)),(TestTrack.bl(1,ind)-orig(1))); 
    else
    c(kk,1) = 1;%-atan2((x((kk-1)*n+3)-orig(2)),(x((kk-1)*n+1)-orig(1)))+atan2((TestTrack.br(2,ind)-orig(2)),(TestTrack.br(1,ind)-orig(1))); % Compute nonlinear inequalities at x.
    end
end

ceq = [];   % Compute nonlinear equalities at x.