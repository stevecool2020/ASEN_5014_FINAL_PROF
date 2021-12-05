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
