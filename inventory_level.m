function [time,demand,level_diff,level] = inventory_level(d,lambda,eta,rho,xi,time0,delta_t,Q)
% generate inventory levels and inventory changes
% input parameter:
% d: basic demand
% lambda, eta: Weibull pars
% rho: sensitivity coefficient of preservation investment
% xi: preservation cost
% time0: the time of order arrival
% delta_t: the time resolution
% Q: the order quantity
% output parameter
% time: sampling time
% demand: demand quantity
% level: inventory level


% calculate order cycle based on order quantity
T=Q/d;
% the moment when the inventory drops to 0
tT=time0+T;
% generate sampling time with time resolution as the step size
time=[(time0+delta_t):delta_t:tT]';
% demand quantity


demand=d*exp(-lambda * exp(-rho*xi)*(time-time0).^eta);
% inventory levels
level=d*exp(-lambda * exp(-rho*xi) * (time-time0).^eta).*(tT-time);
% inventory changes
level_diff=diff([Q;level]);

end

