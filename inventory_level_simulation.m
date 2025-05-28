function [time_simu,demand_simu,level_simu] = inventory_level_simulation(alpha,beta,p,std_dev,lambda0,rho,xi,time0,delta_t,Q)
% simulate inventory levels and inventory changes
% input parameter:
% alpha: basic demand
% beta: price sensitivity coefficient
% p: price
% std_dev: standard deviation of error in demand regression equation
% lambda0: natural deteriorating rate
% rho: parameters in the reduced deterioration rate ratio function
% xi: preservation cost
% time0: the time of order arrival
% delta_t: the time resolution
% Q: the order quantity
% output parameter:
% time_simu: simulated sampling time
% demand_simu: simulated demand quantity
% level_simu: simulated inventory level


% actual deterioration rate
lambda=lambda0*exp(-rho*xi);
% initial inventory data
time_simu=[];
demand_simu=[];
level_simu = [];
% initial inventory level
level_remain=Q;
% record t_{k-1}
time_k1=time0;
% record t_k
time_k=time_k1+delta_t;
% demand quantity
demand = lambda\(alpha-beta*p)*(exp(-lambda.*(time_k1-time0))-exp(-lambda.*(time_k-time0)));
% keep demand >0
demand_error = delta_t*demand + std_dev*randn;
while demand_error<0
    % demand with error = demand + random error
    demand_error = delta_t*demand + std_dev*randn;
end
% deteriorating quantity
deterioration=lambda*(0.5*levelattime(alpha,beta,p,lambda0,rho,xi,time_k,time0,Q) ...
    +0.5*levelattime(alpha,beta,p,lambda0,rho,xi,time_k1,time0,Q));
% reduction expectation
deterioration_poissrnd = poissrnd(delta_t*(deterioration));
while deterioration_poissrnd == 0
    % reduction expectation
    deterioration_poissrnd = poissrnd(delta_t*(deterioration));
end
% random reduction
reduction=deterioration_poissrnd+demand_error;
% level remain = level remain - random reduction
level_remain=level_remain-reduction;
while level_remain > 0
    % store sampling time
    time_simu=[time_simu;time_k];
    demand_simu=[demand_simu;demand_error];
    % store inventory level and changes
    level_simu=[level_simu;level_remain];
    % update t_{k-1}
    time_k1=time_k;
    % update t_{k}
    time_k=time_k+delta_t;
    % demand quantity
    demand = lambda\1*(alpha-beta*p)*(exp(-lambda.*(time_k1-time0))-exp(-lambda.*(time_k-time0)));
    % demand with error = demand + random error
    demand_error = demand + std_dev*randn; % 
    while demand_error<0
        demand_error = demand + std_dev*randn;
    end
    % deteriorating quantity
    deterioration=lambda*(0.5*levelattime(alpha,beta,p,lambda0,rho,xi,time_k,time0,Q) ...
        +0.5*levelattime(alpha,beta,p,lambda0,rho,xi,time_k1,time0,Q));
    % reduction expectation
    % disp(delta_t*(deterioration))
    deterioration_poissrnd = poissrnd(delta_t*(deterioration));
    % while deterioration_poissrnd == 0
    %     % reduction expectation
    %     deterioration_poissrnd = poissrnd(delta_t*(deterioration));
    % end
    % random reduction
    reduction=deterioration_poissrnd+demand_error;
    % level remain = level remain - random reduction
    level_remain=level_remain-reduction;
end


end

