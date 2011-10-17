
function xdot = diff_eq(x,t)


%dipole = .01*exp(-(t-3)^2/2);
%E_in = .1*exp(-t^2);
[kappa, gamma] = constants ();
xdot = zeros (2,1);

xdot(1) = dipole(t)*x(2)-gamma*x(1); %sigma

xdot(2) = -dipole(t)*x(1) - kappa*x(1)+sqrt(2*kappa)*E_in(t);  %E

end
