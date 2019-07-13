function [C, next_pos] = random_positions(init_pos, prec_pos, n, m)
%2 consecutive random positions generation
% Input arguments are: 
% Initial position
% Precedence position
% Dimension of system n 
% Number of agents m

% Initialization
C = zeros(m,n);
final_pos = zeros(m,n);
diff_pos = zeros(m,n);

for i = 1:m
    for j = 1:n
        % I add or subtract 0, 1 or 2 from the initial position to obtain the second position of my agent
        final_pos(i,j) = prec_pos(i,j) + randi([-2,2]); 
      
        % I compute the difference between initial position and final position of my agent
        diff_pos(i,j) = final_pos(i,j)- init_pos(i,j);
        
        % I create the output matrix C = [C1' C2' .. Cm']'
        C(i,j) = -2*diff_pos(i,j);
    end    
end

% figure('Name','Initial position (blue), Final position (red)')
% scatter(target_pos(1,1),target_pos(1,2), 'x', 'g')
% axis([-15 15 -15 15])
% hold on
% 
% for i = 1:m 
%     scatter(init_pos(i,1),init_pos(i,2), 'b')
%     scatter(final_pos(i,1),final_pos(i,2), 'r')    
% end
next_pos = final_pos;

end

