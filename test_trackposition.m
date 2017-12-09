%% Use this file to check the nonlinear track constraint generation
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
flag=1; %1 to test points inside track 0 to test points outside track
if flag==1
    load('Adam_ips.mat','Z')
    xv_mat=Z([1000:1000:8000],1); yv_mat=Z([1000:1000:8000],3);%Points inside track
elseif flag==0
    xv_mat=[230,460,575,860,1008,1235]; yv_mat=[20,230,350,515,490,509];
end

for kk=1:length(xv_mat)
    %     xv=287; yv=-176;
    xv=xv_mat(kk); yv=yv_mat(kk);
    %Find closest points on track
    [mindist,ind]=sort(((xv-TestTrack.cline(1,:)).^2+(yv-TestTrack.cline(2,:)).^2),'ascend');
    %Locations of closest centerline points
    xc1=TestTrack.cline(1,ind(1)); yc1=TestTrack.cline(2,ind(1));
    xc2=TestTrack.cline(1,ind(2)); yc2=TestTrack.cline(2,ind(2));
    mc=(yc2-yc1)/(xc2-xc1);%Slope of centerline segment
    %Locations of closest left track points
    xl1=TestTrack.bl(1,ind(1)); yl1=TestTrack.bl(2,ind(1));
    xl2=TestTrack.bl(1,ind(2)); yl2=TestTrack.bl(2,ind(2));
    ml=(yl2-yl1)/(xl2-xl1);%Slope of left track segment
    cl=yl2-ml*xl2;%Yaxis intersection for left track segment
    %Locations of closest right track points
    xr1=TestTrack.br(1,ind(1)); yr1=TestTrack.br(2,ind(1));
    xr2=TestTrack.br(1,ind(2)); yr2=TestTrack.br(2,ind(2));
    mr=(yr2-yr1)/(xr2-xr1);%Slope of right track segment
    cr=yr2-mr*xr2;%Yaxis intersection for right track segment
    
    %Location of point at closest perpendicular distance from vehicle on
    %centerline
    xc=xc1+(1/(1+mc^2))*((xv-xc1)+mc*(yv-yc1));
    yc=yc1+(mc/(1+mc^2))*((xv-xc1)+mc*(yv-yc1));
    
    %Intermediate variables (bl and br should be >1 and <-1 for car to lie inside
    %track)
    br=(mr*xc+cr-yc)/((yv-yc)-mr*(xv-xc));
    bl=(ml*xc+cl-yc)/((yv-yc)-ml*(xv-xc));
    
    %Intermediate variables (invbl and invbr should be between -1 and 1  for car to lie inside
    %track)
    invbr=((yv-yc)-mr*(xv-xc))/(mr*xc+cr-yc);
    invbl=((yv-yc)-ml*(xv-xc))/(ml*xc+cl-yc);
    
    sign(invbr)*invbr-1
    sign(invbl)*invbl-1
    %Stop if the constraint is wrong (violated when it shouldn't be)
    if flag==1
        if(sign(invbr)*invbr-1)>0
            keyboard
        end
        if(sign(invbl)*invbl-1)>0
            keyboard
        end
    elseif flag==0
        if(sign(invbr)*invbr-1)<0
            keyboard
        end
        if(sign(invbl)*invbl-1)<0
            keyboard
        end
    end
    %Point on left track corresponding to current vehicle position
    xli=xc+bl*(xv-xc);
    yli=yc+bl*(yv-yc);
    
    %Point on right track corresponding to current vehicle position
    xri=xc+br*(xv-xc);
    yri=yc+br*(yv-yc);
    
    figure(h1)
    hold on
    %     plot(xv,yv,'om')
    %     plot(xc,yc,'vm')
    plot([xc xv],[yc yv],'-vm')
    plot([xc xri],[yc yri],'-om')
    plot([xc xli],[yc yli],'-om')
end
