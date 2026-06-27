function [time_simu,demand_simu,level_diff_simu,level_simu] = inventory_level_simulation(d,std_dev,lambda,eta,rho,xi,time0,delta_t,Q)
% simulate inventory levels and inventory changes
% input parameter:
% d: basic demand
% std_dev: standard deviation of error in demand regression equation
% lambda,eta:  Weibull pars
% rho: parameters in the reduced deterioration rate ratio function
% xi: preservation cost
% time0: the time of order arrival
% delta_t: the time resolution
% Q: the order quantity
% output parameter:
% time_simu: simulated sampling time
% demand_simu: simulated demand quantity
% level_diff_simu: simulated inventory change
% level_simu: simulated inventory level


% initial inventory data
time_simu=[];
demand_simu=[];
level_diff_simu = [];
level_simu = [];
% initial inventory level
level_remain=Q;
% record t_{k-1}
time_k1=time0;
% record t_k
time_k=time_k1+delta_t;


% demand quantity
demand = d*(exp(-lambda*exp(-rho*xi).*((time_k1+time_k)/2)^eta));
% disp(demand)
% keep demand >0
demand_error = delta_t*demand + std_dev*randn;
while demand_error<0
    % demand with error = demand + random error
    demand_error = delta_t*demand + std_dev*randn;
end
% actual deterioration rate
theta_t = lambda * eta * ((0.5*time_k+0.5*time_k1)-time0) ^(eta - 1)  *exp(-rho*xi);
% disp(theta_t)
% deteriorating quantity
deterioration=theta_t*levelattime(d,lambda,eta,rho,xi,0.5*time_k+0.5*time_k1,time0,Q);
% reduction expectation
deterioration_poissrnd = poissrnd(delta_t*(deterioration));
while deterioration_poissrnd == 0
    % reduction expectation
    deterioration_poissrnd = poissrnd(delta_t*(deterioration));
end
% disp(deterioration_poissrnd)
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
    level_diff_simu = [level_diff_simu;-reduction];
    % update t_{k-1}
    time_k1=time_k;
    % update t_{k}
    time_k=time_k+delta_t;
    % demand quantity
    demand = d*(exp(-lambda*exp(-rho*xi).*((time_k1+time_k)/2)^eta));
    % disp(demand)
    % demand with error = demand + random error
    demand_error = demand + std_dev*randn; % 
    while demand_error<0
        demand_error = demand + std_dev*randn;
    end
    % actual deterioration rate
    theta_t = lambda * eta * ((0.5*time_k+0.5*time_k1)-time0) ^(eta - 1)  *exp(-rho*xi);
    % disp(theta_t)
    % deteriorating quantity
    % When the actual I (t) is 0, it must be stopped, otherwise det will be 0 later, which does not conform to the model assumption
    deterioration=theta_t*levelattime(d,lambda,eta,rho,xi,0.5*time_k+0.5*time_k1,time0,Q);
    if deterioration==0
        break;
    end
    % reduction expectation
    deterioration_poissrnd = poissrnd(delta_t*(deterioration));
    while deterioration_poissrnd == 0
        % reduction expectation
        deterioration_poissrnd = poissrnd(delta_t*(deterioration));
    end
    % disp(deterioration_poissrnd)
    % random reduction
    reduction=deterioration_poissrnd+demand_error;
    % level remain = level remain - random reduction
    level_remain=level_remain-reduction;
end


end

