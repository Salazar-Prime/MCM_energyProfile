function pbest = fitness(pos_par,n_par,search_dim)

cost_func = [1,2,1,1 ...    %solar
             1,1,1,1 ...    % wind
             1,1,5,1];  %Geothermal
         
% pos_par = zeros(search_dim + 4,n_par);
% pos_par(1:search_dim,1:n_par) = pos_par_main;
% for i=1:n_par
%     for j=1:4
%         pos_par(j+8,i) = 0.8 - pos_par(j+4,i) + pos_par(j,i);
%     end
% end
% 
% for i=1:n_par
%     for j=1:4
%         temp(i) = 0.8 - pos_par(j+4,i) - pos_par(j,i) - pos_par(j+8,i) ;
%     end
% end

for i=1:n_par
        temp(i) = 0.8 - sum(pos_par(:,i));
end

% pbest calculation
pbest = 0.4*(cost_func * pos_par) + abs(0.6*temp);

end


% x = linspace(-10,10,1000);
% for i=1:1000
% for j = 1:1000
% x1 = x(i); x2 = x(j);
% pbest(i,j) = 0.5 + ( (sin(sqrt(x1.*x1+x2.*x2))).^2 - 0.5 ) ./ ( 0.0001 * (x1.*x1+x2.*x2) + 1 ).^2;
% end
% end