clear all
close all

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

sys = ss(A,B,C,D);

eigVal = eig(A);
[eig_Vec, eig_ValDiag] = eig(A);
figure
plot(eigVal, 'o');
xlabel('Real Part')
ylabel('Imaginary Part')

%%% STUB %%%


%%%% SCRIPT STUB %%%%




%%% FUNCTIONS STUB %%%%
