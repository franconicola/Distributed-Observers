clc;
clear;
close all;

n = 2;

A = 0.1*[-5 1; -1 -5];
%A = [1 0.5; -0.5 2];
B = [ 0; 0.787];
C = [1 0];

sys = ss(A,B,C,0);
eig(A)
% System A,B,C (Dynamic) put in the companion form
%sys_dyn = canon(sys,'companion')
sys_dyn = sys;
if charpoly(sys_dyn.A) == minpoly(sys_dyn.A)
    disp('The minimal polynomial is equal to the characteristic polynomial')
end

% Check the observability by the vector g' (Cdyn)
Odyn = obsv(sys_dyn.A,sys_dyn.C);

if n == rank(Odyn)
    disp("System (A, Cdyn) is Observable")
end

% Characteristic polynomial of A
alpha = charpoly(sys_dyn.A);

% Matrices L and phi creation
l = 1;
L = zeros(n+l,1);

for i= 1:(n+l)

    if(i == 1)
        L(i,:) =  sys_dyn.C*sys_dyn.B;
    else
        for j = 1:i

            if j == 1
                L(i,:) = sys_dyn.C*(sys_dyn.A^(i-1))*sys_dyn.B;
            else
                L(i,:) = L(i,:) + alpha(j)*sys_dyn.C*(sys_dyn.A^(i-j))*sys_dyn.B;
            end
        end
    end
end  

phi = [1 alpha(2:n+1);  
       zeros(1,1) L(1:n,:)';
       L(1:n+l,:)']';

rank(phi)

%desired coefficients
beta = charpoly(-10*eye(n+l));
fprintf('Desired polynomial coefficients:\n');
disp(beta);

gamma = beta(1,2:n+1) - alpha(1,2:n+1);
gamma = [gamma beta(1,n+2)];

X = linsolve(phi,gamma')

Bdyn = [sys_dyn.B zeros(n,l); zeros(l,l) eye(l)];
Cdyn = [sys_dyn.C zeros(1,l); zeros(l,n) eye(l)];
Kdyn = [-X(3)  -X(2)+X(3)*X(1); 1 -X(1)];


Ktr = Bdyn*Kdyn;
Abar = Kdyn(2,2);
Bbar = Kdyn(2,1);
Cbar = [ 0; Ktr(2,2)];
Dbar = Ktr(1:2,1);

Hdyn = [sys_dyn.A+Dbar*sys_dyn.C Cbar; Bbar*sys_dyn.C Abar];
Hdyn1 = [sys_dyn.A zeros(2,1); zeros(1,3)] + Bdyn*Kdyn*Cdyn;


charpoly(Hdyn)   
Bdyn = [-Dbar; -Bbar];
sys_dyn = ss(Hdyn,Bdyn,eye(3),zeros(3,1));

pzmap(sys_dyn)


T = 10;
dt = 0.001;
t = 0:dt:T; 

u = 1*(t);
x0 = [1;0];
sys = ss(A,B,eye(2),0);
x = step(sys,t);  %,t,x0);
y = C*x';
x_obs = lsim(sys_dyn,y,t,zeros(3,1));

error = x - x_obs(:,1:2);

figure 
plot(t,x, 'b')
hold on
plot(t,x_obs(:,1:2),'r')
plot(t,x_obs(:,3),'g')
figure
plot(t,error)

