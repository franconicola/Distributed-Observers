clc;
clear;
close all;


%% DEFINITION OF SYSTEM

%Number of agents
m = 3;
%Number of state variables
n = 4;
%Increased dimension of 3° observer
n_ = m-1;

%State Matrix
A = [0 1 -2 0; -1 0 3 0; -2 0 0 2; 0 0 -2 1];

%Input Matrix
B = eye(n);

%Output Matrix
C = eye(n);

%Output Matrices for each observers
C1 = [1 0 0 0; 0 1 0 0];
C2 = [0 1 0 0; 0 0 1 0];
C3 = [0 0 1 0; 0 0 0 1];
CJ = [C1' C2' C3']';

% Conditions for Observers

Q = obsv(A,CJ);

if rank(A) == rank(Q) 
    disp("System (A,C) is jointly Observable")
else
    disp("System (A,C) is not Observable")
end 

sys = ss(A,B,C,0);

x0 = [1.0 5.0 2.0 0.0];
t = 0:0.01:10;

%initial(sys,x0,t)
pole(sys)



%% DEFINITION OF MATRICES FOR OBSERVERS

%Kronecker product for matrix A
A_O = kron(eye(m),A);

%The transpose of incidence matrix of communication graph 
E = [0 1 0; 1 0 1; 0 1 0]';

e1 = E(1,1:3);
e2 = E(2,1:3);

%Kronecker product for vectors e1, e2
C12 = kron(e1,eye(n));
C23 = kron(e2,eye(n));

%Matrices defined by unit vector in position k
b1 = [1; 0; 0];
b2 = [0; 1; 0];
b3 = [0; 0; 1];

%Kronecker product for matrices b1, b2 and b3
B1 = kron(b1, eye(n));
B2 = kron(b2, eye(n));
B3 = kron(b3, eye(n));






%% CALCULATION OF STABILIZABLE MATRICES

poles = [-1+1*1i -1-1*1i -1 -1]*2; %desired poles

R = eye(n)*0.1;

%k = place(A,B,poles)'; %state feedback controller gain matrix
k = lqr(A,B,C,R)'

K1 = k(1:2,1:4)';
K2 = k(2:3,1:4)';
K3 = k(3:4,1:4)';

%state feedback controller gain matrix

obs_poles = [1+1*1i 1-1*1i 1 1]; %desired poles

F = place(A',C',obs_poles)'

F12 = F'; 
F21 = F12;
F32 = F12;
F23 = F32;

%Observer 1:
H11 = -B1*K1*C1*B1';
H12 = B1*F12*C12;

%Observer 2:
H22 = -B2*K2*C2*B2';
H21 = B2*F21*C12;
H23 = B2*F23*C23;

%Observer 3:
H33 = -B3*K3*C3*B3';
H32 = B3*F32*C23;

H = A_O + H11+H12+H22+H21+H23+H33+H32


P1 = ctrb(H,B1);
P2 = ctrb(H,B2);
P3 = ctrb(H,B3);

if rank(H) == rank(P1) && rank(H) == rank(P2) && rank(H) == rank(P3)
    disp("System (H,B1),(H,B2),(H,B3) are Controllables")
else
    disp("System (H,B1),(H,B2),(H,B3) are not Controllables")
end 

Q1 = obsv(H,C12);
Q2 = obsv(H,C23);

if rank(H) == rank(Q1) && rank(H) == rank(Q2)
    disp("System (H,C12) and (H,C23) are Observables")
else
    disp("System (H,C12) and (H,C23) are not Observables")
end 

sys_obs = ss(H,B2,C12,0);

pole(sys_obs)







%%DEFINITION OF OBSERVERS
N = zeros(n);


%Observer 1:
H11 = H(1:4,1:4);
H12 = H(1:4,5:8);
M11 = [1 0 0 0;0 1 0 0; 0 0 0 0;0 0 0 0];
M12 = [0 0 0 0;0 0 0 0; 0 0 1 0;0 0 0 1];

%Observer 2:
H22 = H(5:8,5:8);
H21 = H(5:8,1:4);
H23 = H(5:8,9:12);
M21 = [1 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0];
M22 = [0 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 0];
M23 = [0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 1];

%Observer 3:
H33 = H(9:12,9:12);
H32 = H(9:12,1:4);
M33 = [0 0 0 0;0 0 0 0; 0 0 1 0;0 0 0 1];
M32 = [1 0 0 0;0 1 0 0; 0 0 0 0;0 0 0 0];

%I take in to account the increased dynamic
B_ = -0.2*[1 1; 1 1; 1 1];
D_ = -0.2*ones(n);
C_ = -0.2*[1 1 1; 1 1 1; 1 1 1; 1 1 1];
A_ = -0.2*ones(m);


K3 = [K3; -B_]
H32 = [H32 - D_; -B_ -B_]
H33 = [H33 + D_ C_; B_ B_ A_]
M33 = [M33 [0 0 0; 0 0 0; 0 0 0; 0 0 0]];


sys_obs_3 = ss(H33,[H32 K3*C3],M33, [M32 N]);

pole(sys_obs_3)









