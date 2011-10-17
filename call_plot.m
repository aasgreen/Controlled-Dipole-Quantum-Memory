function [ t,x ,Ef, analy2] = call_plot(amp)

%call_plot.m will solve the coupled ODE's for sigma and E. 
%   The coupled ODE's for sigma and E are in the function dy. That is
%also where you specify gamma and kappa and g_0.

%Main Script to run:

%First, generate the time vector and the initial conditions: 
tspan = [0 5000];
y0 = [0,0]
[kappa, gamma, g_0, c] = constants();
global dip_amp;
dip_amp = amp;
 sig_analy2(4000)
sig_analy(4000)
'Displaying Cavity Mode test'
Elec_analy(800)
Elec_analy2(800)
%Call the solver- ydot has the system of equations.
[t,x] = ode45(@ydot, tspan, y0);

%Start generating the quantities to be plotted.
Eout = zeros( [length(t), 1]);
Einplt = zeros( [length(t),1] );
Dipolplt = zeros ([length(t),1] );
%analy = zeros ([length(t),1]);
analy2 = zeros([length(t), 1]);
%now create plotting data sets
Eout_an = zeros( [length(t), 1]);
Elec_an = zeros( [length(t), 1]);
Elec_an2 = zeros( [length(t), 1]);
for j=1:length(t),
  Eout(j) = -E_in( t(j) ) + sqrt(2*kappa)*x(j,2);
  Einplt(j) = E_in( t(j) );
  Dipolplt(j) = dipole(t(j));
  %analy(j) = sig_analy(t(j));
  analy2(j) = sig_analy2(t(j));
  Eout_an(j) = Elec_out_analy(t(j));
  Elec_an(j) = Elec_analy(t(j));
  Elec_an2(j) = Elec_analy2(t(j));
end 
%plot(t,analy(j))
%Now, create the plots of interest.
%sigma will be stored within x(:,1)
%E will be stored within x(:,2) 
%the (:,1) notation means 'all the rows in colulmn 1'
figure(1)
cla
hold on
%set(1,'defaultaxesfontsize',15)
plot(t,x(:,1),t,imag(x(:,1)) ,t, x(:,2),t,imag(x(:,2)))
h = legend('Re \sigma','Im \sigma','Re E','Im E',1)
hold off
xlabel('ns','fontsize',14)
ylabel('Amplitude','fontsize',14)
%set(h,'fontsize',13)
%print('-dpdf', strcat('sige',num2str(dip_amp),'.pdf'))

figure(2)
cla
hold on
plot(t,Einplt,'--' ,t, Dipolplt,'-.', t, Eout )
h2 = legend('E_{in}', 'N^{1/2}g_0(t) GHz', 'E_{out}',1)
%set(h2,'fontsize',13);
hold off
xlabel('ns', 'fontsize',14)
ylabel('Amplitude', 'fontsize',14)
xlim([0,5000])
%print('-dpdf',strcat('badruninout',num2str(dip_amp),'.pdf'))

figure(3)
plot(t,x(:,1).*conj(x(:,1)),'--', t,x(:,2).*conj(x(:,2)));
h3 = legend('|\sigma|^2', '|E|^2')
xlabel('ns','fontsize',14)
ylabel('Amplitude','fontsize',14)
%print('-dpdf',strcat('badrunexcite',num2str(dip_amp),'.pdf'))

figure(8)
cla
plot(t, imag(analy2), t, imag(x(:,1)));
h10 = legend('Analytical', 'Numerical')
xlabel('ns', 'fontsize',14)
ylabel('Amplitude','fontsize',14)
%Now, we will calculate the efficieny of our memory:

figure(5)
cla
plot(t,Elec_an, t,x(:,2))
h11 = legend('Analytical Cavity (approx)', 'Numerical Cavity')
xlabel('ns', 'fontsize',14)
ylabel('Amplitude','fontsize',14)
cutoff = 1650 %the time that seperates readin from readout
in_index = find(t<cutoff);
out_index = find(t>cutoff);
Ein_total = trapz(t,abs(Einplt).^2)
Ef_in = trapz(t(in_index),abs(Eout(in_index)).^2)/Ein_total
Ef_out =trapz(t(out_index), abs(Eout(out_index)).^2)/Ein_total



 %This function contains the equations to be solved.
    function dy = ydot(t,y)

	%Set the physical constants
	
	[kappa, gamma, g_0, c] = constants();
	g_0 = 1;	
	%initialize the vector to hold the solutions
        dy = zeros(2,1);
	
	%these are the equations to solve
        dy(1) = sqrt(-1)*g_0*dipole(t)*y(2)-gamma*y(1); %sigma

        dy(2) =  -kappa*y(2)+sqrt(2*kappa)*E_in(t) +g_0*sqrt(-1)*dipole(t)*y(1);

    end

end

%function dipole. This will be a chosen function that will be given
%for the dipole. Probably a gaussian
function out = dipole (time)
global dip_amp;
out = dip_amp*exp(-(time-490).^2/((2*210)^2))-dip_amp*exp(-(time-2500).^2/((2*110)^2))-dip_amp*exp(-(time-4000).^2/((2*110)^2));
%out = sin( (time-120)*pi/100);
%out = time^100;
%out = time/100-time/10*sqrt(-1)
end

%function E_in. This is the function for the incoming light
%pulse. Probably a gaussian
function out = E_in (time)
out =  1*exp(-(time-700).^2/(2*200^2));
end

function [kappa, gamma, g_0, c] = constants()
kappa = .9;gamma = 1/200000;
c = 3.0*10^(8);
w_0 = c/765*10^(-9);
h_bar = 1.05457*10^(-34);
V = 10^(-6);
epsilon_0 = 8.854*10^(-12);
g_0 = sqrt(w_0/(2*epsilon_0*h_bar*V));
end



%function sig_analy. This will be the analytical function of sigma. The
%objective will be to test if this matches up with the numerical work.
function out = sig_analy(time)
[kappa, gamma, g_0, c] = constants();
steps = 100
X1 = zeros(steps+1,1);
X2=zeros(steps+1,1);
time = time;
time_array = (0:time/steps:time);
for i = 1:steps+1,
	time_i = time_array(i);
	for j = 1:i,
	time_j = time_array(j);
	  X1(j) = (gamma+dipole(time_j)*dipole(time_j)/kappa);
	end
  X2(i) =sqrt(-1*2/kappa)*dipole(time_i)*E_in(time_i)*exp(trapz(time_array,X1) );
end
RHS =trapz(time_array,X2)
%Now we have calculated everything under the integral,
%we still need to calculate the exp function in the front.
X3 = zeros(100,1);
for i = 1:101,
	X3(i) = -(gamma+dipole(time_array(i))*dipole(time_array(i))/kappa);
end
out = exp(trapz(time_array,X3))*trapz(time_array,X2);
exp(X3);
end

function out = sig_analy2(time);
[kappa, gamma, g_0, c] = constants();
steps = 500; %Number of steps in integration

%calculate the first term within the inner integral
Y1 = [0:time/steps:time;(gamma+(dipole(0:time/steps:time).^2)/kappa)];

%This will create a vector that contains the cumulative
%integration of the values of Y1 at each step
Y2 = cumtrapz(Y1(1,:), Y1(2,:));

%Now, we can caluclate the term that is the under the primary
%integral.
Y3 = sqrt(-2/kappa).*dipole(0:time/steps:time).*E_in(0:time/steps:time).*exp(Y2);
RHS =trapz(Y1(1,:),Y3);%integrate it with respect to time (Y1)

%Now, calculate the lefthand side
Y4 = [0:time/steps:time;-(gamma+(dipole(0:time/steps:time).^2)/kappa)];
out=exp(trapz(Y4(1,:), Y4(2,:)))*RHS;
end

%this function is within the approximation of adiabatic elimantion.
function out = Elec_analy(time);
[kappa, gamma, g_0, c] = constants();
out = sqrt(-1)/kappa*dipole(time).*sig_analy2(time)+sqrt(2/kappa)*E_in(time);
end

function out =  Elec_out_analy(time);
[kappa, gamma, g_0, c] = constants();
out = sqrt(2*kappa)*Elec_analy(time)-E_in(time);
end

%this function solves the E without the adiabatic approx
function out = Elec_analy2(time);
%I don't think we can vecorize, as we are calling sig_analy2, which is
%vecorized, and passing a vector to a vecoriztion won't work well I 
%don't think. Will use loops- be slow, but should work.

[kappa, gamma, g_0, c] = constants();
%number of steps in integration
steps = 800; 
Y1 = [0:time/steps:time; exp(kappa*(0:time/steps:time)).*(sqrt(2*kappa)*E_in(0:time/steps:time)+sqrt(-1)*dipole(0:time/steps:time).*sig_analy2(0:time/steps:time))];

out =1/2*exp(-kappa*time)*trapz(Y1(1,:), Y1(2,:));
end

