function [H,K] = distributed_observer(A,C,E,n,m)
% Distributed observer design function will provides you the matrix H
% Candidates arguments are: 
% plant matrices A and C
% incidence graph matrix E
% dimension of the plant n
% number of agents m

%Kronecker product for matrix A
Abig = kron(eye(m),A);

%The transpose of incidence matrix of communication graph 
Ebig = E';         
Ebig =  kron(Ebig, eye(n));

%Matrices defined by unit vector in position k
Bbig = kron(eye(n), eye(m));

%state feedback controller gain matrix
% poles = -1*ones(1,n);           %desired poles
% K = place(A',C',poles);
K = 0.2*ones(n,n);

% Matrix associate to communication graph
F = kron(eye(n), eye(m));

% Overall observer matrix
H_temp = Abig;

for i = 1:m
       
    if( i < 2)
        k = i;
    else
        k = (i-1)*n +1;
    end
    
    H_temp = H_temp + Bbig(:,k:n*i)*F(k:n*i,k:n*i)*Ebig(k:n*i,:) -Bbig(:,k:n*i)*K(i,:)'*C(i,:)*Bbig(:,k:n*i)';
        
end

H = H_temp;

end

