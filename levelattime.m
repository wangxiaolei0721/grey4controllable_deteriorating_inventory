function level = levelattime(d,lambda,eta,rho,xi,time,time0,Q)
% calculate inventory levels at a specific time
% input parameter:
% d: basic demand
% lambda, eta: Weibull pars
% rho: sensitivity coefficient of preservation investment
% xi: preservation cost
% time0: the time of order arrival
% delta_t: the time resolution
% Q: the order quantity
% output parameter
% level: inventory level


% calculate order cycle based on order quantity
T=Q/d;
% the moment when the inventory drops to 0
tT=time0+T;
% equation
level=d*exp(-lambda * exp(-rho*xi) * (time-time0).^eta).*(tT-time);


end
