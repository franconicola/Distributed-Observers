clear global;
close all;


A = [0 1 0 0; -1 0 0 0; 0 0 0 2; 0 0 -2 0];

B = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

C = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

C1 = C(1:2,1:4);
C2 = C(2:3,1:4);
C3 = C(3:4,1:4);

B1 = B(1,1:4);
B2 = B(2,1:4);
B3 = B(3,1:4);
B4 = B(4,1:4);

P = ctrb(A,B);

if rank(A) == rank(P)
    disp("System is Controllable")
else
    disp("System is not Controllable")
end 

Q = obsv(A,C);

if rank(A) == rank(Q) 
    disp("System is Observable")
else
    disp("System is not Observable")
end 

poles = [-3+3*1i -3-3*1i -3 -3]; %desired poles

K=place(A,B,poles); %state feedback controller gain matrix

K1 = K(1:2,1:4)';
K2 = K(2:3,1:4)';
K3 = K(3:4,1:4)';

%The transpose of incidence matrix of communication graph 
E = [1 0; 1 1; 0 1; 0 0]';

E1 = E(1,1:4);
E2 = E(2,1:4);


%state feedback controller gain matrix

obs_poles = [5 5 5 5]; %desired poles

F = place(A',C',obs_poles)';

F1 = F(1,1:4)';
F2 = F(2,1:4)';
F3 = F(3,1:4)';
F4 = F(4,1:4)';

%Matrix 


%Observer 1:
H11 = A -B1*K1*C1*B1';
H12 = (B1*F1*E1)';
M11 = [1 0 0 0;0 1 0 0; 0 0 0 0;0 0 0 0];
M12 = [0 0;0 0;1 0;0 1]; % 0 0 0 0;0 0 1 0;0 0 0 1];
N1 = [0;0;0;0];

%Observer 2:
H22 = A -B2*K2*C2*B2';
H21 = (B2*F2*E1)';
H23 = (B2*F3*E2)';
M21 = [1; 0; 0; 0]; %0 0 0 0; 0 0 0 0;0 0 0 0];
M22 = [0 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 0];
M23 = [0; 0; 0; 1];%0 0 0 0; 0 0 0 0;0 0 0 1];
N2 = [0 0;0 0;0 0;0 0];

%Observer 3:
H33 = A -B3*K3*C3*B3';
H32 = (B3*F4*E2)';
M32 = [1 0;0 1;0 0;0 0]; %0 1 0 0; 0 0 0 0;0 0 0 0];
M33 = [0 0 0 0; 0 0 0 0;0 0 1 0;0 0 0 1]; 
N3  = [0;0;0;0];


A1 = kron(eye(3),A);
B11 = kron(B1,eye(4));
E11 = kron(E1,eye(4));


R1 = [1 0 0 0];
R2 = [0 1 0 0]; 
R3 = [0 0 1 0];
R4 = [0 0 0 1];


%  t = 0:0.01:2;
%  u = zeros(size(t));
%  x0 = [0.01 0 0];
%  
%  
%  sys = ss(A,B,C,0);
%  
%  
%  
%  [y,t,x] = lsim(sys,u,t,x0);
%  plot(t,y)
%  title('Open-Loop Response to Non-Zero Initial Condition')
%  xlabel('Time (sec)')
%  ylabel('Ball Position (m)')

