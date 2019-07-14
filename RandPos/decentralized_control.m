function [Hdyn,B_obs] = decentralized_control(H,C,E,Kgain,n,m)
% DECENTRALIZED CONTROL
% Input arguments are: 
% Distributed observer matrix H
% Plant matrix C
% Graph matrix E
% Gain matrix Kgain defined in distributed observer H
% Dimension of system n 
% Number of agents m

% Overall dimension of the system
v = n*m;

% Observability and Controllability of the total system from last agent
agent = m;

% Define the size of matrix C for one agent
mes = size(C,1)/m;

% Dynamic Compensator of dimension l implemented to the last observer
l = m-1;

%Matrices defined by unit vector in position k
Bbig = kron(eye(n), eye(m));

%The transpose of incidence matrix of communication graph         
Ebig =  kron(E', eye(n));

%Iterator
k = (agent-1)*n +1;

%Matrix C defined as the matrix between connections
Cm = Ebig(agent,:);     

%We check observability and controllability from the last agent
Obar = obsv(H,C(1+(agent-1)*mes:mes+(agent-1)*mes,:)*Bbig(:,k:n*agent)'); % Cm
Rbar = ctrb(H,Bbig(:,k:n*agent));

if rank(H) <= rank(Obar) && rank(H) <= rank(Rbar)
    disp("System (Cm, H, Bm) is Observable and Controllable")

    if size(charpoly(H)) == size(minpoly(H))
        if charpoly(H) == minpoly(H)
            disp('The minimal polynomial is equal to the characteristic polynomial')
        end
    else
        disp('The minimal polynomial is not equal to the characteristic polynomial')
    end

    % Vector g' used to obtain the matrix phi 
    g = C(mes+(agent-1)*mes,:)*Bbig(:,k:n*agent)';      %Cm

    % Check the observability by the vector g' (Cdyn)
    Odyn = obsv(H,g);

    if v == rank(Odyn)
         disp("System (H, g) is Observable")
    else
         disp("System (H, g) is not Observable")
    end

    % Characteristic polynomial of H
    alpha = charpoly(H);

    % System definition %
    csys_dyn = ss(H,Bbig(:,k:n*agent),C(1+(agent-1)*mes:mes+(agent-1)*mes,:)*Bbig(:,k:n*agent)',zeros(mes,n));

    % Matrices L and phi creation
    L = zeros(v,n);
    phi = zeros(v+l,v+l);

    for i= 1:v
        if(i == 1)
            L(i,:) =  g*csys_dyn.B;
        else
            for j = 1:i+1-l

                if j == 1
                    L(i,:) = g*(csys_dyn.A^(i-1))*csys_dyn.B;
                else
                    L(i,:) = L(i,:) + alpha(j)*g*(csys_dyn.A^(i-j))*csys_dyn.B;
                end
            end
        end
    end  
    
    for i= 1:l+m
        if(i == 1)
            phi(i,:) =  [1 alpha(2:v+1) zeros(1,l-i)];
        elseif(i < l)
            phi(i,:) =  [zeros(1,l-i) 1 alpha(2:v+1) zeros(1,l-i)];
        elseif(i == l)
            phi(i,:) =  [zeros(1,i-1) 1 alpha(2:v+1)];
        elseif(i == l+1)
            phi(l+1:l+n,:) = [zeros(n, l) L'];  
        elseif(i > l+1 && i < m+l)
            phi((l+1+(i-l-1)*n):(l+n+(i-l-1)*n),:) = [zeros(n, 2*l+1-i) L' zeros(n, i-l-1)];
        else
            phi((l+1+(i-l-1)*n):(l+n+(i-l-1)*n),:) = [L' zeros(n, i-l-1)];
        end     
    end

    phi = phi';

    % desired coefficients 
    beta = charpoly(-2*eye(v+l));

    % subtraction between desired coefficients and given by H
    gamma = beta(1,2:v+1) - alpha(1,2:v+1);
    gamma = [gamma beta(1,v+2:v+l+1)];
    
    % solutions delta of the system phi*delta=gamma
    X = linsolve(phi,gamma')

    % Gain implemented in the last observer given by the delta solutions
    Kdyn = zeros(n+l, mes+l);

    for i=1:n+l
        for j=1:mes+l
            if(i <= n) 
                if(j == 1)
                    Kdyn(i,j) = -X(l+n*l+i);
                elseif(j <= mes && n < mes+l)
                    Kdyn(i,j) = 0;
                elseif(j < mes+l)
                    Kdyn(i,j) = -X(l+(j-l)*n+i)+X(l+n*l+i)*X(j-l);
                else
                    Kdyn(i,j) = -X(l+i)+X(l+n*l+i)*X(j-l);
                end
            elseif(i == n+1)
                if(j == 1)
                    Kdyn(i,j) = 1;
                elseif(j <= mes && n < mes+l)
                    Kdyn(i,j) = 0;
                else
                    Kdyn(i,j) = -X(j-mes);
                end
            else
                if(j == 1)
                    Kdyn(i,j) = 0;
                elseif(j <= l && n < mes+l)
                    Kdyn(i,j) = 0;
                elseif( j == l+i-n-1)
                    Kdyn(i,j) = 1;
                else
                    Kdyn(i,j) = 0;
                end
            end
        end
    end
    Kdyn
    
     % Increasing dimensionmes of our last observer
    Bdyn = [csys_dyn.B zeros(v,l); zeros(l,n) eye(l)];

    % Definitions of matrices for increased observer
    Ktr  = Bdyn*Kdyn;
    Abar = Kdyn(n+1:n+l,1:l);
    Bbar = Kdyn(n+1:n+l,l+1:l+mes);
    Cbar = Ktr(1:v,mes+1:mes+l);
    Dbar = Ktr(1:v,1:mes);

    % New H matrix with the dynamic compensator
    Hdyn = [csys_dyn.A+Dbar*csys_dyn.C Cbar;  Bbar*csys_dyn.C Abar];
    
    % The input matrix B_obs
    B_obs = zeros(v+l,mes*m);
    for j = 1:m
        if(j==1)
            B_obs(:,1+(j-1)*mes:mes+(j-1)*mes) = [ Kgain(1+(j-1)*n:n+(j-1)*n,1+(j-1)*mes:mes+(j-1)*mes); zeros(n*(m-1)+l,mes)];
        elseif(j == m)
            B_obs(:,1+(j-1)*mes:mes+(j-1)*mes) = [ [zeros(n*(m-1),mes); Kgain(1+(j-1)*n:n+(j-1)*n,1+(j-1)*mes:mes+(j-1)*mes)]-Dbar; -Bbar];
        else
            B_obs(:,1+(j-1)*mes:mes+(j-1)*mes) = [ zeros(n*(j-1),mes); Kgain(1+(j-1)*n:n+(j-1)*n,1+(j-1)*mes:mes+(j-1)*mes); zeros(n*(m-j)+l,mes)];
        end
    end

else
    Hdyn = H;
    
    % The input matrix B_obs
    B_obs = zeros(v,mes*m);
    for j = 1:m
        if(j==1)
            B_obs(:,1+(j-1)*mes:mes+(j-1)*mes) = [ Kgain(1+(j-1)*n:n+(j-1)*n,1+(j-1)*mes:mes+(j-1)*mes); zeros(n*(m-1),mes)];
        elseif(j == m)
            B_obs(:,1+(j-1)*mes:mes+(j-1)*mes) = [ zeros(n*(m-1),mes); Kgain(1+(j-1)*n:n+(j-1)*n,1+(j-1)*mes:mes+(j-1)*mes)];
        else
            B_obs(:,1+(j-1)*mes:mes+(j-1)*mes) = [ zeros(n*(j-1),mes); Kgain(1+(j-1)*n:n+(j-1)*n,1+(j-1)*mes:mes+(j-1)*mes); zeros(n*(m-j),mes)];
        end
    end
    
    disp("System (Cm, H, Bm) is nor Controllable neither Observable")
end
    
end

