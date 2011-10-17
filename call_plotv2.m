function [ t,x ] = call_plot()

clear
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tspan = [0 300];
y0 = [0,.001]

[t,x] = ode45(@ydot, tspan, y0);
%other actions

disp([t,x])

Eout = zeros( [length(t), 1]);
Einplt = zeros( [length(t),1] );
Dipolplt = zeros ([length(t),1] );

%now create plotting data sets
for j=1:length(t),
  %Eout(j) = -E_in( t(j) ) + sqrt(2*kappa)*x(j,2);
  Einplt(j) = x(j,2);
  Dipolplt(j) = dipole(t(j));
end 


figure(1)
cla
hold on
plot(t,x(:,1), 'r', t, x(:,2))
title('Sigma, E_in')
h = legend('Sigma','E_in',1)
hold off
figure(2)
cla
hold on
plot(t,x(:,2), t, Dipolplt, 'r')
title('Ein, Dipole')
h2 = legend('Ein', 'Dipole',1)
hold off
    function dy = ydot(t,y)
        gamma =.03; kappa =40;
        dy = zeros(2,1);
        dy(1) = dipole(t)*y(2)/sqrt(2*kappa)-gamma*y(1); %sigma

        dy(2) =  sqrt(2*kappa)*(-kappa*y(2)/sqrt(2*kappa)+sqrt(2*kappa)*y(2) +dipole(t)*y(1));

    end

end

%function dipole. This will be a chosen function that will be given
%for the dipole. Probably a gaussian
function out = dipole (time)
out = 2.6*exp(-(time-140)^2/((2*12)^2));
%out = sin( (time-120)*pi/100);
%out = time^100;

end

%function E_in. This is the function for the incoming light
%pulse. Probably a gaussian
function out = E_in (time)
out =  1*exp(-(time-150)^2/100);
end
