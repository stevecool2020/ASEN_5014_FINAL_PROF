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
figure
plot(eigVal, 'o');
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

% isControllable = rank(ctrb(OLsys)) == numel(A(:,1)) % yes
% isObservable = rank(obsv(OLsys)) == numel(A(:,1)) % yes

% Names for plots
OLsys.StateName = {'x (radial)';'y (along-track)';'z (cross-track)';...
    'xdot (radial velocity)';'ydot (along-track velocity)';'zdot (cross-track velocity)'};
OLsys.StateUnit = {'meters';'meters';'meters';...
    'meters/seconds';'meters/seconds';'meters/seconds'};
OLsys.OutputName = OLsys.StateName(1:3);
OLsys.OutputUnit = OLsys.StateUnit(1:3);
OLsys.InputName = {'xddot';'yddot';'zddot'};
OLsys.InputUnit = {'meters/seconds^2';'meters/seconds^2';'meters/seconds^2'};

% Open loop poles
% olpoles = pole(OLsys);
% pzmap(OLsys)

% impulse(OLsys) % this makes sense

% Simulate the plant response to the desired initial conditions (no
% disturbance)

% X0 = [0 10000 0 0 0 0];
% X0 = [10000 0 0 0 0 0];
% X0 = [0 0 10000 0 0 0];
% X0 = [10 375 0 0 0.00009 0];
X0 = [0 400 0 0 0.0235 0];
% Target radius
rTgt_m = 6778E3;
% Target mean motion
nTgt = sqrt(muEarth_m3ps2/rTgt_m^3);
% Target orbit period
TTgt_s = 2*pi/nTgt; 

tvec_s = 0:0.1:2*TTgt_s;

% figure();
% initial(OLsys,X0,tvec_s);
[y,t,x] = initial(OLsys,X0,tvec_s); %<----------------this can be ref trajectory
step(OLsys,TTgt_s)

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
despoles = [-0.09 -0.06 -0.03 -0.02 -0.02 -0.01]; %dominant pole is real and at -2 
%% Use place command to get feedback gain:
K = place(A,B,despoles);

% Define closed-loop dynamics
F1 = eye(2); %unity feedforward gain on ref inputs 
F2 = inv(C/(-A+B*K)*B);
F = F2;
Acl = A - B*K;
Bcl = B*F;
Ccl = C;
Dcl = D;
CLsys = ss(Acl,Bcl,Ccl,Dcl);


umax_mps2 = 1E-4; % thrust_kgmps2/massChaser_kg;

rhistvec =  x(:,1:3);

% Check that closed-loop system specs met; change despoles otherwise
%%get response to first reference input profile: 
[Y_CL1,~,X_CL1] = lsim(CLsys,rhistvec,tvec_s);

%compute resulting actuator efforts in each case, where u = -Kx + Fr
U_CL1 = -K*X_CL1' + F*rhistvec'; 

figure("Name",'Manual Pole Placement Actuator Effort')
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

figure("Name","Manual Pole Placement: Response Compared To Desired Position State")
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
plot(tvec_s,rhistvec(:,3),'DisplayName','Reference')
legend show


%% #5 Luenberger Observer 
% todo remove these gains once part 4 is complete
K = zeros(3,6);
F = eye(3);

% check for controllability of observer system
Co_observer = ctrb(A', C');
is_observer_controllable = rank(Co_observer) == nstates

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

figure();
initial(LOCLsys,2*zeros(2*nstates,1),10) %initial observer error states zero
title('Observer error transient responses | Zero intial error')

figure();
initial(LOCLsys,2*ones(2*nstates,1),10) %initial observer error states non-zero
title('Observer error transient responses | With inital error')



%%%% SCRIPT STUB %%%%




%%% FUNCTIONS STUB %%%%
