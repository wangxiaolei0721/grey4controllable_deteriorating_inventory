function profit_average = profit(K,c,h,p,time0,pars,T,xi) 
% profit function
% input parameter:
% K: ordering cost per cycle
% c: production cost
% h: holding cost per unit per unit of time
% p: selling price
% time0: the time of order arrival
% pars: [lambda,eta,rho,d]
% T: replenishment cycle
% xi: preservation cost
% output parameter:
% profit: total profits at T

%
lambda = pars(1);
eta = pars(2);
rho = pars(3);
d0 = pars(4);
%
% Phi = @(T,xi) integral(@(t) exp(-lambda*exp(-rho*xi).*(t-time0).^eta),time0,time0+T);
% Omega = @(T,xi) integral(@(t) exp(-exp(-rho*xi).*(t-time0).^eta).*(time0+T-t),time0,time0+T);


Phi = @(T,xi) arrayfun(@(t,x) ...
    integral(@(s) exp(-lambda*exp(-rho*x).*(s-time0).^eta), ...
    time0, time0+t), T, xi);

Omega = @(T,xi) arrayfun(@(t,x) ...
    integral(@(s) exp(-exp(-rho*x).*(s-time0).^eta).*(time0+t-s), ...
    time0, time0+t), T, xi);

profit_average = ...
    (p*d0.*Phi(T,xi) - K - c*d0.*T - (h+xi).*d0.*Omega(T,xi))./T;

end
