clc;
clear;
close all;

%% SYSTEM'S DEFINITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of agents
m = 3;
%Number of state variables
n = 4;
%Overall dimension of the system         
v = m*n;

%State Matrix
A = [0 1  0 0; 
    -1 0  0 0; 
     0 0  0 2; 
     0 0 -2 0]

%Output Matrices for each observers
C1 = [1 0 0 0]; 
C2 = [0 1 1 0]; 
C3 = [0 0 1 0]; 
      
C = [C1' C2' C3']'

% Associate graph matrix
E = [0 1 0; 
     0 0 1;
     1 0 0];

% Conditions for Observervability
Q = obsv(A,C);

if rank(A) == rank(Q) 
    disp("System (A,C) is jointly Observable")
else
    disp("System (A,C) is not Observable")
end 

% Check the strongly connected components
GraphE = digraph(E);
%plot(GraphE)
s = conncomp(GraphE);

if s == ones(1,m)
    disp("Graph is strongly connected")
else
    disp("Graph is weakly connected")
end

%% DEFINITION OF MATRIX H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Kronecker product for matrix A
Abig = kron(eye(m),A);

%The transpose of incidence matrix of communication graph 
Ebig = E';         
Ebig =  kron(Ebig, eye(n));

%Matrices defined by unit vector in position k
Bbig = kron(eye(n), eye(m));

%state feedback controller gain matrix
poles = -1*ones(1,n); %desired poles          
Kgain = place(A,eye(n),poles);
%Kgain = eye(n);

F = eye(n*m);
K = kron(eye(m), diag(Kgain));

Kbig = [K -F];

%Definition of Cbig
Cbig = zeros(m,v);
for i = 1:m
    
    zeroVectorForward = zeros(1,n*(m-i));
    zeroVectorBackward = zeros(1,n*(i-1));
    if( i > 1 ) 
        Cbig(i,:) = [zeroVectorBackward C(i,1:n) zeroVectorForward];
    elseif (i == 1)
        Cbig(i,:) = [C(i,1:n) zeroVectorForward];
    elseif (i == m)
        Cbig(i,:) = [zeroVectorBackward C(i,1:n)]; 
    end
end

Cbig = [Cbig ; Ebig];

% H matrix
H = Abig -Kbig*Cbig

% Check the controllability and observability of the overall system
% by the last observer
Obar = obsv(H,Cbig(m,1:v));
Rbar = ctrb(H,Bbig(1:v,(m-1)*n+1:v));

if v == rank(Obar) && v == rank(Rbar)
    disp("System (Cbig, H, Bbig) is Observable and Controllable")
    rho_cntr = rank(Bbig(1:v,1:n));
else
    disp("System (Ebig, H, Bbig) is nor Controllable neither Observable")
    rho_cntr = 0;
end

%% DECENTRALIZED CONTROL THEORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dynamic Compensator of dimension l implemented to the last observer
l = m-1;

% System creates for the last observer 
csys_dyn = ss(H,Bbig(1:v,(m-1)*n+1:v),Cbig(m,1:v),zeros(1,n));

% Check if it's cyclic 
if charpoly(csys_dyn.A) == minpoly(csys_dyn.A)
    disp('The minimal polynomial is equal to the characteristic polynomial')
% if not put in the companion form
else
    csys_dyn = canon(csys_dyn,'companion');
end

% Characteristic polynomial of H
alpha = charpoly(csys_dyn.A);

% Matrix L creation
L = zeros(v,n);
for i= 1:v
    if(i == 1)
        L(i,1:n) =  csys_dyn.C*csys_dyn.B;
    else
        for j = 1:i
            if j == 1
                L(i,1:n) = csys_dyn.C*(csys_dyn.A^(i-1))*csys_dyn.B;
            else
                L(i,1:n) = L(i,1:n) + alpha(j)*csys_dyn.C*(csys_dyn.A^(i-j))*csys_dyn.B;
            end
        end
    end
end   

% Matrix phi creation
phi = zeros(v+l,rho_cntr*(l+1) + l);
index = 0;
for i = 1:l+rho_cntr-1
    if i <= l
        phi(i,:) = [zeros(1,i-1) 1 alpha(2:v+1) zeros(1,l-i)];
    else       
        phi(i+index*(rho_cntr-1):i+(index+1)*(rho_cntr-1),:) = [zeros(rho_cntr,2*l+1-i) L' zeros(rho_cntr,i-l-1)];
        index = index + 1;
    end
    
end
phi = phi';
% Desired polynomial coefficients 
beta = charpoly(-1*eye(v+l));
fprintf('Desired polynomial coefficients:\n');
disp(beta)

% Subtraction between beta and alpha
gamma = beta(1,2:v+1) - alpha(1,2:v+1);
gamma = [gamma beta(1,v+l:v+l+1)]';

% Delta coefficients by solution of phi*X=gamma system
X = linsolve(phi,gamma)

% Increasing dimension of system's matricies
Bdyn = [csys_dyn.B zeros(v,l); zeros(l,n) eye(l)];
Cdyn = [csys_dyn.C zeros(1,l); zeros(l,v) eye(l)];

% K gain given by the compensator
Kdyn = zeros(rho_cntr+l,1+l);
index = 2;
for i=1:rho_cntr+l
    if i <= rho_cntr
        Kdyn(i,:) = [-X(i+l+l*rho_cntr) -X(i+l+(l-1)*rho_cntr) -X(i+l)+X(i+l+l*rho_cntr)*X(2)+X(i+l+(l-1)*rho_cntr)*X(1)]; 
    elseif i == rho_cntr + 1
        Kdyn(i,:) = [1 0 -X(2)];
    else
        Kdyn(i,:) = [0 1 -X(1)];
    end
end
Kdyn
% Stabilization of the overall system
Hdyn = [csys_dyn.A zeros(v,l);  zeros(l,v+l)] + Bdyn*Kdyn*Cdyn;
        
fprintf('Obtained polynomial coefficients:\n');
disp(charpoly(Hdyn))

% % Increasing dimension of system's matricies
Bbig_dyn = [Bbig zeros(v,l); zeros(l,v) eye(l)];
Cbig_dyn = [Cbig(1:m,:) zeros(m,l); zeros(l,v) eye(l)];
 
sys_dyn = ss(Hdyn,Bbig_dyn,Cbig_dyn,zeros(l+m,v+l));
pzmap(sys_dyn) 


T = 10;
dt = 0.001;
t = 0:dt:T; 

%u = ones(4,1)*(t);
x0 = 0.1*ones(14,1);
sys = ss(A,ones(4,1),eye(4),zeros(4,1));
x = step(sys,t);
y = C*x'; 
x_obs = lsim(sys_dyn,[y; y; y; y; zeros(2,4)*x'],t,zeros(14,1));

error = x - x_obs(:,1:4);

figure 
plot(t,x, 'b')
hold on
plot(t,x_obs(:,1:4),'r')
plot(t,x_obs(:,5),'g')
figure
plot(t,error)


%% DEFINITION OF OBSERVERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Observer 1:
H11 = H(1:n,1:n);
H13 = H(1:n,2*n+1:v);
K1 = K(1:n,1);
M11 = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
M13 = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

N = zeros(n);

%Observer 2:
H22 = H(n+1:2*n,n+1:2*n);
H21 = H(n+1:2*n,1:n);
K2 = K(n+1:2*n,2);
M22 = [0 0 0 0; 0 1 0 0;0 0 1 0;0 0 0 0];
M21 = [1 0 0 0;0 0 0 0; 0 0 0 0;0 0 0 1];

%Observer 3:
H33 = Hdyn(2*n+1:v+l,2*n+1:v+l);
H32 = Hdyn(2*n+1:v+l,n+1:2*n);
K3 = [K(2*n+1:v,3); zeros(l,1)];
M33 = [0 0 0 0 0 0;0 0 0 0 0 0; 0 0 1 0 0 0;0 0 0 1 0 0];
M32 = [1 0 0 0;0 1 0 0; 0 0 0 0;0 0 0 0];

















