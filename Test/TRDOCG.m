%% TRDOCG 
% Test for Random generation of Distributed Observers with different 
% Graph Connectivity 

clc;
clear;
close all;

%% DEFINITION OF SYSTEM

%Number of agents
m = 4;
%Number of state variables
n = 6;
%Overall dimension of the system         
v = m*n;

%Number of state variables seen by each observer
numberOfOnes_C = int8(n/4);
%Maximum number of connections in the graph for each observer
numberOfOnes_E = int8(m/2);

% Get a list of random locations, with no number repeating.
for i = 1:m
    indexes_C = randperm(n);
    indexes_E = randperm(m);
    
    % Start off with all zeros.
    c = zeros(1, n);
    e = zeros(1, m);
    
    % Now make half of them, in random locations, a 1.
    c(indexes_C(1:numberOfOnes_C)) = 1;
    e(indexes_E(1:numberOfOnes_E)) = 1;
    
    %Output Matrix
    C(i,:) = c;
    
    %Graph Matrix
    E(i,:) = e; % ones(1,m);                    
    E(i,i) = 0;
end

%State Matrices
A = randi([0, 1], n)        
C

%% Check the Jointly Observability
O = obsv(A,C);

if rank(A) == rank(O) 
    disp("System is Jointly Observable")
else
    disp("System is not Jointly Observable")
end 

%% Check the strongly connected components
GraphE = digraph(E);
plot(GraphE)
s = conncomp(GraphE);

if s == ones(1,m)
    disp("Graph is strongly connected")
else
    disp("Graph is weakly connected")
end


%% DEFINITION OF MATRICES FOR OBSERVERS

%Kronecker product for matrix A
Abig = kron(eye(m),A);

%The transpose of incidence matrix of communication graph 
Ebig = E';         
Ebig =  kron(Ebig, eye(n));

%Matrices defined by unit vector in position k
Bbig = kron(eye(n), eye(m));

%state feedback controller gain matrix
poles = -2*ones(1,n);                   %desired poles
Kgain = place(A',eye(n),poles);         %only in the principal diagonal

F = eye(n*m);
K = kron(eye(m), diag(Kgain));

Kbig = [K -F];

%Definition of C_
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

H = Abig -Kbig*Cbig;

%% CHECK
% I am checking if the H creates by the previous steps, is equal than the H
% defined in the paper
H_check = Abig;

%Observers and Neighbours
for i = 1:m
       
    if( i < 2)
        k = i;
    else
        k = (i-1)*n +1;
    end
    
    S = Bbig(1:v,k:n*i)*F(k:n*i,k:n*i)*Ebig(k:n*i,1:v) -Bbig(1:v,k:n*i)*K(k:i*n,i)*C(i,1:n)*Bbig(1:v,k:n*i)';
        
    H_check = H_check + S;
end

H_check = H_check - H;

if H_check == zeros(v,v)
    disp("H finds with H = Abig -Kbig*Cbig it is equal than the H defined in the paper")
end

%% Observability of the total system
Obar = obsv(Abig,Cbig);

if rank(Abig) == rank(Obar) 
    disp("System (Abig,Cbig) is Observable")
    
else
    disp("System (Abig,Cbig) is not Observable")
end 


%obs_poles = -2*ones(1,v); %desired poles
%Kgain_obs = place(Abig',Cbig',obs_poles);

%H = Abig - Kgain_obs*Cbig;


%Observability form
Bbar = eye(v);         %kron(eye(m,m),ones(n,1));
Dbar = zeros(v+m,v);

%[Abar,Bbar,Cbar,T,k] = obsvf(Abig,Bbar,Cbig);
sys = ss(H,Bbar,Cbig,Dbar);
%pzmap(sys)

