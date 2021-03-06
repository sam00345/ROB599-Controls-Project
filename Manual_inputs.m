%Manually define inputs that complete the track

load('TestTrack.mat')

bl_x = TestTrack.bl(1,:);
bl_y = TestTrack.bl(2,:);

br_x = TestTrack.br(1,:);
br_y = TestTrack.br(2,:);

cline_x = TestTrack.cline(1,:);
cline_y = TestTrack.cline(2,:);

theta = TestTrack.theta(1,:);

U1 = [-0.005*ones(900,1) 1600*ones(900,1)];
U2 = [-0.19*ones(665,1) 0*ones(665,1)];
U3 = [0*ones(800,1) 0*ones(800,1)];
U4 = [0.06*ones(375,1) 1000*ones(375,1)];
U5 = [-0.2*ones(500,1) 600*ones(500,1)];
U6 = [0.1*ones(425,1) 0*ones(425,1)];
U7 = [0*ones(300,1) 800*ones(300,1)];
U8 = [-0.045*ones(825,1) 300*ones(825,1)];
U9 = [-0.001*ones(375,1) 0*ones(375,1)];
U10 = [0*ones(25,1) -9000*ones(25,1)];
U11 = [0.01*ones(300,1) 0*ones(300,1)];
U12 = [0.15*ones(400,1) 0*ones(400,1)];
U13 = [-0.010*ones(700,1) 200*ones(700,1)];
U14 = [-0.12*ones(600,1) 300*ones(600,1)];
U15 = [0.1*ones(400,1) 0*ones(400,1)];
U16 = [0.0012*ones(650,1) 1750*ones(650,1)];
U17 = [0.08*ones(250,1) 0*ones(250,1)];
U18 = [0.18*ones(135,1) 0*ones(135,1)];
U19 = [0*ones(800,1) 4000*ones(800,1)];

U = [U1;U2;U3;U4;U5;U6;U7;U8;U9;U10;U11;U12;U13;U14;U15;U16;U17;U18;U19];

Z = forwardIntegrateControlInput(U); 

plot(bl_x,bl_y,br_x,br_y,cline_x,cline_y,'--',Z(:,1),Z(:,3),'k')