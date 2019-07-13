clc;
clear;
close all;


%% DEFINITION OF SYSTEM

%Number of agents
m = 3;
%Number of state variables
n = 4;

%State Matrix
A = [0 1 0 0; -1 0 0 0; 0 0 0 2; 0 0 -2 0];

%Input Matrix
B = eye(n);

%Output Matrix
C = eye(n);

%Output Matrices for each observers
C1 = [1 0 0 0; 0 1 0 0];
C2 = [0 1 0 0; 0 0 1 0];
C3 = [0 0 1 0; 0 0 0 1];
CJ = [C1' C2' C3']';

%Check the Jointly Observability
O = obsv(A,C);

if rank(O) == n 
    disp("System is Jointly Observable")
else
    disp("System is not Jointly Observable")
end 

sys = ss(A,B,C,0);
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
B1 = kron(b1, eye(n))
B2 = kron(b2, eye(n));
B3 = kron(b3, eye(n));

%% CALCULATION OF STABILIZABLE MATRICES

poles = [-1 -1 -1 -1]*20; %desired poles

k = place(A,B,poles)' %state feedback controller gain matrix

K1 = k(1:2,1:4)';
K2 = k(2:3,1:4)';
K3 = k(3:4,1:4)';

%state feedback controller gain matrix

obs_poles = [1 1 1 1]; %desired poles

F = place(A',CJ',obs_poles)';

F12 = F(1:4,1:4)';
F21 = F12;
F32 = F(1:4,3:6)';
F23 = F32;

%Observer 1:
H11 = -B1*K1*C1*B1'
H12 = B1*F12*C12

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

%%DEFINITION OF OBSERVERS

%Observer 1:
H11 = H(1:n,1:n);
H12 = H(1:4,5:8);
M11 = [1 0 0 0;0 1 0 0; 0 0 0 0;0 0 0 0];
M12 = [0 0 0 0;0 0 0 0; 0 0 1 0;0 0 0 1];

N = zeros(n);

%Observer 2:
H22 = H(5:8,5:8);
H21 = H(5:8,1:4);
H23 = H(5:8,9:12);
M21 = [1 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 0];
M22 = [0 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 0];
M23 = [0 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 1];

%Observer 3:
H33 = H(9:12,9:12);
H32 = H(9:12,5:8);
M33 = [0 0 0 0;0 0 0 0; 0 0 1 0;0 0 0 1];
M32 = [1 0 0 0;0 1 0 0; 0 0 0 0;0 0 0 0];










