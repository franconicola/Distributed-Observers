function [H,K] = distributed_observer(A,C,E,n,m)
% Distributed observer design function will provides you the matrix H
% Candidates arguments are: 
% plant matrices A and C
% incidence graph matrix E
% dimension of the plant n
% number of agents m

%Define the size of matrix C for one agent
mes = size(C,1)/m;

%Kronecker product for matrix A
Abig = kron(eye(m),A);

%The transpose of incidence matrix of communication graph          
Ebig =  kron(E', eye(n));

%Matrices defined by unit vector in position k
Bbig = kron(eye(n), eye(m));

%state feedback controller gain matrix
K = zeros(n*m,m*mes);

for i = 1:m
    K(1+(i-1)*n:n+(i-1)*n,1+(i-1)*mes:mes+(i-1)*mes) = 0.01*randi([1,9],[n,mes]);
end


% Matrix associate to communication graph
F = kron(eye(n), eye(m));

% Overall observer matrix
H_temp = Abig;

for i = 1:m
    %iterator   
    k = (i-1)*n +1;
 
    H_temp = H_temp + Bbig(:,k:n*i)*F(k:n*i,k:n*i)*Ebig(k:n*i,:) -Bbig(:,k:n*i)*K(1+(i-1)*n:n+(i-1)*n,1+(i-1)*mes:mes+(i-1)*mes)*C(1+(i-1)*mes:mes+(i-1)*mes,:)*Bbig(:,k:n*i)';
        
end

H = H_temp;

end

