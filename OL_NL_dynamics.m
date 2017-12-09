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
%% Forward integrate nonlinear dynamics with input 
x0 = [287 5 -176 0 2 0]';
u{1}=[0*ones(300,1) 5000*ones(300,1)];%
u{2}=[-0.015*ones(50,1) 3000*ones(50,1)];%
u{3}=[-0.02*ones(50,1) 2000*ones(50,1)];%
u{4}=[-0.03*ones(50,1) 4000*ones(50,1)];%
u{5}=[-0.025*ones(100,1) -1000*ones(100,1)];%
u{6}=[-0.032*ones(100,1) -1000*ones(100,1)];%
u{7}=[ -0.05*ones(100,1) -1000*ones(100,1)];%
u{8}=[-0.07*ones(100,1) -1000*ones(100,1)];%
u{9}=[ -0.15*ones(100,1)  -1000*ones(100,1)];%
u{10}=[-0.3*ones(300,1)  -2000*ones(300,1)];%
u{11}= [0*ones(50,1)  -500*ones(50,1)];
u{12}=[0.1*ones(50,1)  500*ones(50,1)];
u{13}=[0*ones(50,1)  1000*ones(50,1)];
u{14}=[-0.05*ones(150,1)  0*ones(150,1)];
u{15}=[0.2*ones(100,1)  0*ones(100,1)];
u{16}=[-0.05*ones(100,1)  2000*ones(100,1)];
u{17}=[-0.001*ones(350,1)  3000*ones(350,1)];
u{18}=[0.00*ones(100,1)  -5000*ones(100,1)];
u{19}=[0.11*ones(100,1)  -2000*ones(100,1)];
u{20}=[0.2*ones(100,1)  0*ones(100,1)];
u{21}=[0*ones(250,1)  0*ones(250,1)];
u{22}=[0.1*ones(100,1)  50*ones(100,1)];
% u{23}=[0.*ones(100,1)  500*ones(100,1)];

% u{23}=[0.2*ones(100,1)  500*ones(100,1)];
% % u{23}=[0.0*ones(50,1)  1000*ones(50,1)];

u_all=[];
dt=0.01;
x=x0;
for pp=1:length(u)    
u_all=[u_all; u{pp}];
% if pp>1
%     ucurr=[u{pp-1}(end,:);u{pp}(1:end-1,:)];
% else
%     ucurr=[u{pp}(1:end-1,:)];
% end
ucurr=u{pp};
[Y]=nldynamics(ucurr,x',dt);
Y=Y';
x=Y(:,end);


figure(h1)
hold on
plot(Y(1,:),Y(3,:),'Color',rand(1,3),'Linewidth',1);   
end

[Yall]=nldynamics(u_all,x0',dt);
Yall=Yall';
figure(h1)
hold on
plot(Yall(1,:),Yall(3,:),'-.m','Linewidth',2);   
Yall(:,end)