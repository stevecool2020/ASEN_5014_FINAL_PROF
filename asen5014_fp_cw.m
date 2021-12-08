% asen5014_fp.m
clc; close all;

%% Inputs, Parameters, and Constants
% Earth mean radius:
rEarth_m = 6371E3;
muEarth_m3ps2 = 398600*(1000^3);

% Target radius
rTgt_m = 6778E3;
% Target mean motion
nTgt = sqrt(muEarth_m3ps2/rTgt_m^3);
% Target orbit period
TTgt_s = 2*pi/nTgt; 

% Acceleration equations (linearized)
xdd = @(yd,x) 2*nTgt*yd + 3*nTgt^2*x;
ydd = @(xd) -2*nTgt*xd;
zdd = @(z) -nTgt^2*z;

% Linearized State Space Model
A = [zeros(3),eye(3)
   3*nTgt^2 0 0 0 2*nTgt 0
   0 0 0 -2*nTgt 0 0
   0 0 -nTgt^2 0 0 0];
B = [zeros(3);eye(3)];
C = [eye(3),zeros(3)];
D = zeros(3);

%% Open loop analysis
olsys = ss(A,B,C,D);



% Names for plots
olsys.StateName = {'x (radial)';'y (along-track)';'z (cross-track)';...
    'xdot (radial velocity)';'ydot (along-track velocity)';'zdot (cross-track velocity)'};
olsys.StateUnit = {'meters';'meters';'meters';...
    'meters/seconds';'meters/seconds';'meters/seconds'};
olsys.OutputName = olsys.StateName(1:3);
olsys.OutputUnit = olsys.StateUnit(1:3);
olsys.InputName = {'xddot';'yddot';'zddot'};
olsys.InputUnit = {'meters/seconds^2';'meters/seconds^2';'meters/seconds^2'};

% Open loop poles
% olpoles = pole(olsys);
% pzmap(olsys)

% Simulate the plant response to the desired initial conditions (no
% disturbance)
% step(olsys)
% impulse(olsys) % this makes sense


% X0 = [0 10000 0 0 0 0];
% X0 = [10000 0 0 0 0 0];
% X0 = [0 0 10000 0 0 0];
% X0 = [10 375 0 0 0.00009 0];
X0 = [0 400 0 0 0.0235 0];
tvec_s = 0:0.1:2*TTgt_s;
% figure();
% initial(olsys,X0,TTgt_s);
[y,t,x] = initial(olsys,X0,tvec_s);
step(olsys,tvec_s)

figure();
subplot(311)
plot(t,y(:,1),'DisplayName','Radial Obs');  hold on;
plot(t,x(:,1),'DisplayName','Radial State'); 
legend('show'); ylabel('meters'); grid minor;
subplot(312)
plot(t,y(:,2),'DisplayName','Along-Track Obsv'); hold on;
plot(t,x(:,2),'DisplayName','Along-Track State'); 
legend('show'); ylabel('meters'); grid minor;
subplot(313)
plot(t,y(:,3),'DisplayName','Cross-Track Obsv'); hold on;
plot(t,y(:,3),'DisplayName','Cross-Track State');
legend('show'); ylabel('meters'); grid minor;
xlabel('Time (seconds)');
sgtitle('Initial Conditions System Response')

%% Infinite Horizon Controller


% define umax 
thrust_kgmps2 = 25;
massChaser_kg = 100;

umax_mps2 = thrust_kgmps2/massChaser_kg;
rmag = 0.2;
% rhistvec1 = [rmag.*ones(size(tvec_s));zeros([2,numel(tvec_s)])]';
% rhistvec2 = [zeros([1,numel(tvec_s)]);rmag.*ones(size(tvec_s));zeros([1,numel(tvec_s)])]';
% rhistvec3 = [zeros([2,numel(tvec_s)]);rmag.*ones(size(tvec_s))]';
rhistvec1 = [x(:,1)';zeros([2,numel(tvec_s)])]';
rhistvec2 = [zeros([1,numel(tvec_s)]);x(:,2)';zeros([1,numel(tvec_s)])]';
rhistvec3 = [zeros([2,numel(tvec_s)]);x(:,3)']';

awts = ones([1,numel(A(:,1))]);
rho = 1;
awts = awts./sum(awts);


%Q = [eye(3)*10,zeros(3);zeros(3),eye(3)*1];% eye(size(A))*10;
%R=eye(numel(B(1,:)));
Q = diag(awts./[100 100 100 0.01 0.01 0.01].^2);
R = rho*diag(1./(umax_mps2.*ones([1,3])).^2);

[Ks, W, E] = lqr(A,B,Q,R);
% usLqr = nan([3, numel(x(:,1))]);
% 
% for i = 1:numel(x(:,1))
%    usLqr(:,i) = -R^-1 * B'*W*x(i,:)';
% end
% 
% figure()
% subplot(311)
% plot(t,x(:,4),'DisplayName','Radial State'); hold on
% plot(t,usLqr(1,:),'r','DisplayName','Radial InfHorCtrl'); 
% legend('show'); ylabel('m/s'); grid minor;
% subplot(312)
% plot(t,x(:,5),'DisplayName','Along-Track State'); hold on
% plot(t,usLqr(2,:),'r','DisplayName','Along-Track InfHorCtrl'); 
% legend('show'); ylabel('m/s'); grid minor;
% subplot(313)
% plot(t,x(:,6),'DisplayName','Cross-Track State'); hold on
% plot(t,usLqr(3,:),'r','DisplayName','Cross-Track InfHorCtrl'); 
% legend('show'); ylabel('m/s'); grid minor;
% xlabel('Time (seconds)');
% sgtitle('Actuator Effort for Desired Conditions Tracking')

F = (C/(-A+B*Ks)*B)^-1;
Acl = A - B*K;
Bcl = B*F;
Ccl = C;
Dcl = D;
CLsys = ss(Acl,Bcl,Ccl,Dcl);

[ycl1,~,xcl1] = lsim(CLsys,rhistvec1,tvec_s,[X0(1:3),0,0,0]);
[ycl2,~,xcl2] = lsim(CLsys,rhistvec2,tvec_s,[X0(1:3),0,0,0]);
[ycl3,~,xcl3] = lsim(CLsys,rhistvec3,tvec_s,[X0(1:3),0,0,0]);
ucl1 = -K*xcl1' + F*rhistvec1';
ucl2 = -K*xcl2' + F*rhistvec2';
ucl3 = -K*xcl3' + F*rhistvec3';

figure("Name",'Actuator Effort')
subplot(311), hold on
plot(tvec_s,ucl1(1,:))
plot(tvec_s,ucl1(2,:))
plot(tvec_s,ucl1(3,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
subplot(312), hold on
plot(tvec_s,ucl2(1,:))
plot(tvec_s,ucl2(2,:))
plot(tvec_s,ucl2(3,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
subplot(313), hold on
plot(tvec_s,ucl3(1,:))
plot(tvec_s,ucl3(2,:))
plot(tvec_s,ucl3(3,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');

figure("Name","Response Compared To Desired Position State")
subplot(311)
plot(tvec_s,ycl1(:,1),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec1(:,1),'DisplayName','Reference');
legend show
subplot(312)
plot(tvec_s,ycl2(:,2),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec1(:,2),'DisplayName','Reference');
legend show
subplot(313)
plot(tvec_s,ycl1(:,3),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec1(:,3),'DisplayName','Reference')
legend show


%% 





