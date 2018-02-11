%% Varun Aggarwal - 14110142

%% PSO optimization (main function)

clear all; close all
% initializing parameters
% rng(45312)4
n_par = 500000;
search_dim = 12;
c1 = 0.5; c2 = 0.5;
max_iter = 500;
wmax = 1; wmin = 0.8;
% func = @(x1,x2) 0.5 + ((sin(sqrt(x1.^2+x2.^2))).^2-0.5)./((0.001.*(x1.^2+x2.^2)+1).^2); 

% initializing particles
max_pos = 0.2;
max_vel = 0.001;
pos_par = Randomizer(0,max_pos,search_dim,n_par);
vel_par = Randomizer(-max_vel,max_vel,search_dim,n_par);

% inital values of pbest and gbest
pbest = fitness(pos_par,n_par,search_dim);
pbest_coordinates = pos_par;
gbest = min(pbest);
gbest_coordinates = pos_par(:,find(pbest==gbest));

f = figure;
% iterations
for i=1:max_iter
    w = wmax - (i/max_iter)*(wmax-wmin);
    i
    % updating particle
    for j=1:n_par
        % update velocity
%         vel_par(:,j) = w*vel_par(:,j) + c1*rand()*(pbest_coordinates(:,j)-pos_par(:,j)) + c2*rand()*(gbest_coordinates - pos_par(:,j)); 
        temp_v = w*vel_par(:,j) + c1*rand()*(pbest_coordinates(:,j)-pos_par(:,j)) + c2*rand()*(gbest_coordinates - pos_par(:,j)); 
        
        % boundary condition
        for k=1:search_dim
            if abs(temp_v(k)) < max_vel
                vel_par(k,j) = temp_v(k);
            end
        end
        
        % updating position
        temp_p = pos_par(:,j) + vel_par(:,j);
        
        % boundary condition
        for k=1:search_dim
            if temp_p(k) < max_pos && temp_p(k) > 0
                pos_par(k,j) = temp_p(k);
            end
        end
    end % end update particles
    
    % checking fitness
    pbest_temp = fitness(pos_par,n_par,search_dim);
    
    % updating pbest for each particle
    for j=1:n_par
        if(pbest(j)>pbest_temp(j))
            pbest(j) = pbest_temp(j);
            pbest_coordinates(:,j) = pos_par(:,j);
        end
    end
    
    % updating gbest
    gbest_temp = min(pbest);
    if gbest>gbest_temp
        gbest = gbest_temp(1);
        temp = find(pbest==gbest);
        gbest_coordinates = pbest_coordinates(:,temp(1));
    end
    
%        % plotting
%        x = linspace(-10,10,100);
%        for k = 1:100
%            for j = 1:100
%                schaffer(k,j) = func(x(k),x(j));
%            end
%        end

    final_gbest(i) = gbest;
%     hold off
%     contour(x,x,schaffer)
%     hold on
%     plot(gbest_coordinates(1),gbest_coordinates(2),'*r')
%     title(strcat('Gbest over Schaffer function:',num2str(i),'   gbest:',num2str(gbest)))
%     xlabel('x1');ylabel('x2')
%     pause(0.1)
end

plot(final_gbest)
% comment to toggle saving video file
% v = VideoWriter('Trial_6.avi');
% open(v);
% writeVideo(v,MM);
% close(v);  

