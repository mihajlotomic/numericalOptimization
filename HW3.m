% 
% 
x_k = [0;-1];
g = Grad_Fx(x_k(1), x_k(2));
fx = Fx(x_k(1),x_k(2));
B = Hess_Fx(x_k(1),x_k(2));
% B is positive def for x = (0,-1)
%       _     _ 
% B  = | 42 0  |   
%      | 0  20 | 
%       -      - 

% ( B+lambda * eye(2) ) * p_star = -g
% When lambda = 0 => B*p_star = -g; 
% Delta bound by ||p|| = inv(B) x -g
% When lambda > 0; need to find lambda

% Define a range for delta within the defined region
% Then build our p_star data
delta = .25:0.25:2;
p_star = zeros(2,length(delta));

if x_k(2) == -1
    count = 0;
    lambda_0 = 1;
    for delta = .25:0.25:2
        if delta > 0 && delta < norm(B\-g)
            lambda = newton_lambda(delta, B, g, lambda_0); %lambda variable this region
            count = count + 1;
            p_star(:,count) = [2/(42+lambda); 20/(20+lambda)];
        elseif delta >= norm(B\-g) && delta <= 2
            count = count + 1;
            p_star(:,count) = B\(-g);
        else
            error('Invalid delta');
        end
    end
else % x_k = [0, 0.5]
    count = 1;
    lambda_0 = 18.5;
    for delta = .25:0.25:2
        lambda = newton_lambda(delta, B, g, lambda_0); %lambda variable this region
        p_star(:,count) = [2/(-18+lambda); -10/(20+lambda)];
        count = count + 1;
    end
end

%plot the contours
%first build the grids for x1,x2
x1 = -1:.05:1;
x2 = -1:.05:1;
[X1,X2] = meshgrid(x1,x2);
%the planes created by meshgrid are not X1,X2
%feed these in to find the contours
%we know: Fx, g, Hess_Fx - what about p?
%p is just a single representation on the grid
%p_star is a 2 x length(delta)
%what is the size of M and how do we create a p of the right size?

FX =  Fx(X1,X2);
% p = repmat(p_star,
GT = Grad_Fx(X1,X2)';
delta_count = 1;
for delta_slice = 1:length(p_star)
    for i = 1:length(X1)
        count = 1;
        for j = 1:2:length(X2)*2
            GTP(i,count,delta_count) = GT(i,j:j+1)*p_star(:,delta_slice);
            count = count +1;  
        end    
    end
    delta_count = delta_count + 1;
end

% - 402x402 matrix per the size of our meshgrid
B = Hess_Fx(X1,X2);
delta_count = 1;

for delta_slice = 1:length(p_star)
    row_count = 1;
    for i = 1:2:length(X1)*2
        column_count = 1;
        for j = 1:2:length(X2)*2
            PTBP(row_count,column_count,delta_count) = p_star(:,delta_slice)'*B(i:i+1,j:j+1)*p_star(:,delta_slice);
            column_count = column_count +1;  
        end  
        row_count = row_count + 1;
    end
    delta_count = delta_count + 1;
end
FX_D = repmat(FX, [1 1 size(GTP,3)]); % Create a 201x201x9 - 
M = FX_D + GTP + 0.5.*PTBP;
for i = 1:size(M,3)
    meshc(X1,X2,M(:,:,i))
    hold on;
end

% Now plot the p_star values
p_star_slice = p_star(:,1:end);
for i = 1:length(p_star_slice)
    %plot([0, p_star_slice(1,i)], [0,p_star_slice(2,i)],'->');
    quiver(0,0,p_star_slice(1,i), p_star_slice(2,i),'r');
    hold on;
end
