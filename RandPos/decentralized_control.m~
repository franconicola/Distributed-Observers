function [Hdyn,Bbar,Dbar] = decentralized_control(H,C,n,m)
% DECENTRALIZED CONTROL
% Input arguments are: 
% Distributed observer matrix H
% Plant matrix C
% Dimension of system n 
% Number of agents m

% Overall dimension of the system
v = n*m;

% Observability and Controllability of the total system from last agent
agent = m;

%Matrices defined by unit vector in position k
Bbig = kron(eye(n), eye(m));

if( agent < 2)
    k = agent;
else
    k = (agent-1)*n +1;
end
 
Obar = obsv(H,C(agent,:)*Bbig(:,k:n*agent)');
Rbar = ctrb(H,Bbig(:,k:n*agent));

if v == rank(Obar) && v == rank(Rbar)
    disp("System (Cm, H, Bm) is Observable and Controllable")
else
    disp("System (Cm, H, Bm) is nor Controllable neither Observable")
end

% Dynamic Compensator of dimension l implemented to the last observer
l = m-1;

% Vector g' used to obtain the matrix phi 
g = C(agent,:)*Bbig(:,k:n*agent)';


if charpoly(H) == minpoly(H)
    disp('The minimal polynomial is equal to the characteristic polynomial')
end

% Check the observability by the vector g' (Cdyn)
Odyn = obsv(H,g);

if v == rank(Odyn)
     disp("System (H, g) is Observable")
else
     disp("System (H, g) is not Observable")
end

% Characteristic polynomial of H
alpha = charpoly(H);

% System definition
csys_dyn = ss(H,Bbig(:,k:n*agent),C(agent,:)*Bbig(:,k:n*agent)',zeros(1,n));

% Matrices L and phi creation
L = zeros(v+l,l+1);
phi = zeros(v+l,(m*(l+1) + l));

for i= 1:(v+l)
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

for i= 1:(m+l)
    if(i == 1)
        phi(i,:) =  [1 alpha(2:v+1) zeros(1,l-i)];
    elseif(i < l)
        phi(i,:) =  [zeros(1,l-i) 1 alpha(2:v+1) zeros(1,l-i)];
    elseif(i == l)
        phi(i,:) =  [zeros(1,l-i) 1 alpha(2:v+1)];
    elseif(i == l+1)
        phi(i:i+m-1,:) = [zeros(m, l) L(1:v,:)'];
    elseif(i > l+1 && i <m+l)
        phi((i-1+(i-l-1)*m):(i-m+(i-l)*m),:) = [zeros(m, 2*l-i) L(1:v,:)' zeros(m, i-l)];
    else
        phi((i-1+(i-l-1)*m):(i-m+(i-l)*m),:) = [L(1:v,:)' zeros(m, i-l-1)];
    end
end

phi = phi';
rank(phi);

% desired coefficients 
beta = charpoly(-2*eye(v+l));
% disp('Desired eigenvalues:');
% disp(eig(-2*eye(v+l)));

% subtraction between desired coefficients and given by H
gamma = beta(1,2:v+1) - alpha(1,2:v+1);
gamma = [gamma beta(1,v+l+1)];

% solutions delta of the system phi*delta=gamma
X = linsolve(phi,gamma');

% Increasing dimension of our last observer
Bdyn = [csys_dyn.B zeros(v,l); zeros(l,n) eye(l)];
Cdyn = [csys_dyn.C zeros(1,l); zeros(l,v) eye(l)];

% Gain implemented in the last observer given by the delta solutions
Kdyn = [-X(4)  -X(2)+X(4)*X(1);
        -X(5)  -X(3)+X(5)*X(1);
        1            -X(1)];

% Definitions of matrices for increased observer
Ktr = Bdyn*Kdyn;
Abar = Kdyn(3,2);
Bbar = Kdyn(3,1);
Cbar = Ktr(1:v,2);
Dbar = Ktr(1:v,1);

% New H matrix with the dynamic compensator
Hdyn = [csys_dyn.A+Dbar*csys_dyn.C Cbar;  Bbar*csys_dyn.C Abar];
    
end

