clc;
clear;
close all;

%% SYSTEM'S DEFINITION

%Number of agents
m = 2;
%Number of state variables
n = 2;
%Overall dimension of the system         
v = m*n;
%Dimension of the dynamic compensator
l = m-1;
%Number of positions iteraction
STEPS = 3;


%State Matrices
A = zeros(n,n);
B = eye(n);

E = ones(m,m)-m*eye(m);                             % Incidence graph matrix

%Check the strongly connected components
GraphE = digraph(E);
%plot(GraphE)
s = conncomp(GraphE);

if s == ones(1,m)
    disp("Graph is strongly connected")
else
    disp("Graph is weakly connected")
end

%% POSITIONS 

% Random position generation
init_pos = randi([-10,10], m,n);
target_pos = randi([-10,10],1,n);
positions = zeros(STEPS*m,n);

% Random position generation 
C_temp = zeros(m,n);
C_tot = zeros(STEPS*m,n*m);

for i = 1:STEPS
    
    % Updating positions of drones and the output matrix C_temp
    if (i == 1)
        [C_temp, positions(i:m,:)] = random_positions(init_pos, init_pos, n, m);
    else
        [C_temp, positions((1+(i-1)*m):(m+(i-1)*m),:)] = random_positions(init_pos, positions(1+(i-1)*m:m+(i-1)*m,:), n, m);
    end
    
    % Store the value of matrix C_temp 
    C_tot(i,:) = reshape(C_temp',1,[]);
    
    % Create C matrix
    C = zeros(i*m,n);
    for j = 1:m
       C((1+(j-1)*i):(i+(j-1)*i),:) = C_tot(1:i,1+(j-1)*n:n+(j-1)*n);
    end
    C
    %Check the Jointly Observability
    O = obsv(A,C);

    if n == rank(O) 
        disp("System is Jointly Observable")

        % Definition of matrices for observers

        [H,Kgain] = distributed_observer(A,C,E,n,m);

        [Hdyn,B_obs] = decentralized_control(H,C,E,Kgain,n,m);  
        disp(eig(Hdyn))
        obs_dyn = ss(Hdyn,B_obs,eye(size(Hdyn,1)),zeros(size(Hdyn,1),i*m));
                      
        % SIMULATIONs
        T = 15;
        dt = 0.01;
        t = 0:dt:T; 

        plant = ss(A,B,eye(n),zeros(n,n));
        u = zeros(n,1)*t;
        x0 = target_pos';
        x = lsim(plant,u,t,x0);

        y = C*x';    %C_temp    

        z = lsim(obs_dyn, y,t,zeros(size(Hdyn,1),1));

        for j = 1:m
            error(:,1+(j-1)*n:n+(j-1)*n) = x - z(:,1+(j-1)*n:n+(j-1)*n);
        end

        figure('Name','First Observer(red) and plant(black)')
        plot(t,x,'k')
        hold on
        plot(t,z(:,1:2), 'r')

        figure('Name','Second Observer(red), plant(black) and compensator dynamic(green)')
        plot(t,x,'k')
        hold on
        plot(t,z(:,3:4), 'r')
        plot(t,z(:,5), 'g')

        figure('Name', 'Errors of observers')
        plot(t, error)
        grid on
        
    else
        disp("System is not Jointly Observable")
    end
    
    pause 

       
end



fprintf('Target Position: X: %i, Y: %i \n', x(1,1), x(1,2));










