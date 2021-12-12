clear; close all

%%%%%% A,B,C,D Matrices %%%%%%
n = sqrt(398600 / 6778^3);
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*n^2 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -n^2 0 0 0];

B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];

C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];

D = [0 0 0;
     0 0 0;
     0 0 0];

nstates = size(A,1); % how many states? 
OLsys = ss(A,B,C,D);

eigVal = eig(A);
[eig_Vec, eig_ValDiag] = eig(A);
figure('name', 'Pole-Zero Map');
pzmap(OLsys); 
xlabel('Real Part')
ylabel('Imaginary Part')


% Earth mean radius:
rEarth_m = 6371E3;
muEarth_m3ps2 = 398600*(1000^3);


%% #2
% Define control system objectives, requirements. Determine open loop plant
% poles. Simulate the plant respone (no controller) to the desired initial
% conditions, reference inputs and /or exogeneous disturbances, as
% appropriate to control objectives, verify that thees responses make sense

isControllable = rank(ctrb(OLsys)) == numel(A(:,1)) % yes
isObservable = rank(obsv(OLsys)) == numel(A(:,1)) % yes

% Names for plots
OLsys.StateName = {'x (radial)';'y (along-track)';'z (cross-track)';...
    'xdot (radial velocity)';'ydot (along-track velocity)';'zdot (cross-track velocity)'};
OLsys.StateUnit = {'meters';'meters';'meters';...
    'meters/seconds';'meters/seconds';'meters/seconds'};
OLsys.OutputName = OLsys.StateName(1:3);
OLsys.OutputUnit = OLsys.StateUnit(1:3);
OLsys.InputName = {'xddot';'yddot';'zddot'};
OLsys.InputUnit = {'meters/seconds^2';'meters/seconds^2';'meters/seconds^2'};


% define umax 
thrust_kgmps2 = 25;
massChaser_kg = 100; 

% Open loop poles
% olpoles = pole(OLsys);
% pzmap(OLsys)

% impulse(OLsys) % this makes sense

% Simulate the plant response to the desired initial conditions (no
% disturbance)

% Examine different initial conditions
% X0 = [10000 0 0 0 0 0]; % radial offset
% X0 = [0 10000 0 0 0 0]; % in-track offset
% X0 = [0 0 10000 0 0 0]; % cross-track offset
% X0 = [10 375 0 0 0.00009 0];
X0 = [100 5000 -1000 0 0 0]; % TODO : finalize these numbers
% Target radius
rTgt_m = 6778E3;
% Target mean motion
nTgt = sqrt(muEarth_m3ps2/rTgt_m^3);
% Target orbit period
TTgt_s = 2*pi/nTgt; 

tvec_s = 0:0.1:2*TTgt_s;

[y,t,x] = initial(OLsys,X0,tvec_s); 

figure('Name','Open Loop Response to Initial Conditions');
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


%% #4 Manual Pole Placement

% Specify desired closed-loop poles
despoles = [-0.05 -0.09 -0.06 -0.03 -0.02 -0.02 -0.01 -0.005 - 0.0025 -0.0001]./10; %dominant pole is real and at -2 
%% Use place command to get feedback gain:
% K = place(A,B,despoles);

Aaug = [A, zeros(6,3);
        -C, zeros(size(C,1))];
Baug = [B; zeros(size(C,1), size(B,2))];
Caug = [C, zeros(3,3)]; 
Daug = zeros(3,3);
Faug = [zeros(size(B));
        eye(3)]; % F using final value theorem
    
Kaug = place(Aaug,Baug,despoles); 
Aaugcl = Aaug - Baug*Kaug;
Baugcl = Faug;
Caugcl = Caug;
Daugcl = Daug;
CLaugsys = ss(Aaugcl,Baugcl,Caugcl,Daugcl);

umax_mps2 = thrust_kgmps2/massChaser_kg;

rhistvec = zeros(3,numel(tvec_s))';

X0aug = [X0, 0, 0, 0]; % set the augmenmted integrated states to zero initial conditions
% Check that closed-loop system specs met; change despoles otherwise
%%get response to first reference input profile: 
[Y_CL1,~,X_CL1] = lsim(CLaugsys,rhistvec,tvec_s, X0aug);

%compute resulting actuator efforts in each case, where u = -Kx + Fr
% U_CL1 = -Kaug*X_CL1' + Faug*rhistvec'; 
U_CL1 = -Kaug*X_CL1'; 


figure('Name','Manual Pole Placement Actuator Effort')
subplot(311), hold on
plot(tvec_s,U_CL1(1,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
subplot(312), hold on
plot(tvec_s,U_CL1(2,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
subplot(313), hold on
plot(tvec_s,U_CL1(3,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');

figure('Name','Manual Pole Placement: Response Compared To Desired Position State')
subplot(311)
plot(tvec_s,Y_CL1(:,1),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec(:,1),'DisplayName','Reference');
legend show
subplot(312)
plot(tvec_s,Y_CL1(:,2),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec(:,2),'DisplayName','Reference');
legend show
subplot(313)
plot(tvec_s,Y_CL1(:,3),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec(:,3),'DisplayName','Reference');
legend show


%% #5 Luenberger Observer 
% todo remove these gains once part 4 is complete
K = zeros(3,6);
F = eye(3);

% check for controllability of observer system
Co_observer = ctrb(A', C');
is_observer_controllable = rank(Co_observer) == nstates;

% New augmented state for luenberger observer: Xaug_obs = [x; ex; ey; ez; exdot; eydot; ezdot]

% Specify desired closed-loop poles for the Luenberger observer
desobsvpoles =  [-12 -11 -10 -9, -8, -7]; % TODO : Fine tune 
% Use place command to get Luenberger observer gain
L = (place(A',C',desobsvpoles))';

% define new closed loop observer system
AaugLOCL = [A B*K;
           zeros(nstates) A-L*C];
% BaugLOOL = [B; B]; % open loop version
BaugLOCL = [B*F; zeros(size(B))];
CaugLOCL = [zeros(nstates,nstates), eye(nstates)]; %define errors as outputs (should --> 0)
DaugLOCL = zeros(6,3);
LOCLsys = ss(AaugLOCL, BaugLOCL, CaugLOCL, DaugLOCL);

figure('Name','Observer Error Transient 0 IC');
initial(LOCLsys,2*zeros(2*nstates,1),10) %initial observer error states zero
title('Observer error transient responses | Zero Intial Error')

figure('Name','Observer Error Transient With IC');
initial(LOCLsys,2*ones(2*nstates,1),10) %initial observer error states non-zero
title('Observer error transient responses | With Inital Error')


%% # 6 Infinite Horizon Controller

umax_mps2 = thrust_kgmps2/massChaser_kg;
rmag = 0.2;
% rhistvec1 = x(:,1:3);
rhistvec1 = zeros(3,numel(tvec_s))';

% Weights 
awts = ones([1,numel(A(:,1))]);
rho  = 10;
awts = awts./sum(awts);

poswts = [20 90 60];
velwts = 1*[1 1 1];
Q = diag(awts./[poswts, velwts].^2);
R = rho*diag(1./(umax_mps2.*[1 1 1]).^2);

[Ks, W, E] = lqr(A,B,Q,R);
% usLqr = nan([3, numel(x(:,1))]);
% 
% for i = 1:numel(x(:,1))
%    usLqr(:,i) = -R^-1 * B'*W*x(i,:)';
% end

Ff = (C/(-A+B*Ks)*B)^-1;
Acl = A - B*Ks;
Bcl = B*Ff;
Ccl = C;
Dcl = D;
CLsys = ss(Acl,Bcl,Ccl,Dcl);

[ycl1,~,xcl1] = lsim(CLsys,rhistvec1',tvec_s,X0);

ucl1 = -Ks*xcl1' + Ff*rhistvec1';

figure("Name",'LQR Actuator Effort')
subplot(311), hold on
plot(tvec_s,ucl1(1,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
ylabel({'Acceleration';'Radial (m/s^2)'});
subplot(312), hold on
plot(tvec_s,ucl1(2,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
ylabel({'Acceleration';'In-Track (m/s^2)'});
subplot(313), hold on
plot(tvec_s,ucl1(3,:))
plot(tvec_s,umax_mps2*ones(size(tvec_s)),'--k');
plot(tvec_s,-umax_mps2*ones(size(tvec_s)),'--k');
ylabel({'Acceleration';'Cross-Track (m/s^2)'});
xlabel('Time (s)')

figure("Name","LQR Response Compared To Desired Position State")
subplot(311)
plot(tvec_s,ycl1(:,1),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec1(:,1),'DisplayName','Reference');
ylabel({'Distance';'Radial (m)'});
legend show
subplot(312)
plot(tvec_s,ycl1(:,2),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec1(:,2),'DisplayName','Reference');
ylabel({'Distance';'In-track (m)'});
legend show
subplot(313)
plot(tvec_s,ycl1(:,3),'DisplayName','Observed'); hold on;
plot(tvec_s,rhistvec1(:,3),'DisplayName','Reference')
ylabel({'Distance';'Cross-track (m)'}); 
xlabel('Time (s)')
legend show



%%%% SCRIPT STUB %%%%




%%% FUNCTIONS STUB %%%%
