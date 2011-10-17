function [ t,x ,Ef,Ef_analy2,a] = call_plot_new(amp,width)

%call_plot_new.m will solve the coupled ODE's for sigma and E. 
%   The coupled ODE's for sigma and E are in the function dy. That is
%also where you specify gamma and kappa and g_0.
%This is the current expiremental program. The main optimized one is
%call_plot.m

%This one will test the analytical predictions to obtain perfect retrieval.
%I will be starting will retrevial.
%Main Script to run:

%First, generate the time vector and the initial conditions: 
tspan = [0 5000];
y0 = [0,0] %retrieval, so |sigma(t=0)| = 1
[kappa, gamma, g_0, c] = constants();
global initial_cond;
initial_cond=y0;
global dip_amp;
dip_amp = amp;
global dip_width;
dip_width = width;
 sig_analy2(4000);
'Displaying Cavity Mode test'
%Elec_analy(800)
%Elec_analy2(800)
%Call the solver- ydot has the system of equations.
[t_out,x_out] = ode45(@ydot, tspan, y0);
 t = t_out';
x = x_out;
%To be consistant with the rest of the program, take transpose
%to get row vectors of the above quantities
Tot_Eout = 0;
Tot_Ein = 0;
%Now we want to calculate the readout efficiency.
%Hopefully, as we normalized the sigma, this should just be the
%the absolute integral of Eout squared.
Eout = sqrt(2*kappa)*x(:,2)'-E_in(t);
Eout2 = E_in(t).*(-1+2*exp(cumtrapz(t,-dipole(t).^2/kappa)));
%figure(100)
%cla
%plot(t,Eout, t, Eout2)
Ef = 1-trapz(t,abs(Eout).^2)/trapz(t',abs(E_in(t)).^2)
Ef_analy = analy_efout(t);
Ef_analy2= analy_ef2out(t);
Ein_tot=trapz(t,abs(E_in(t)).^2);

%Ef_analy
%Ef
%Ef_analy2
a = trapz(t,dipole(t).^2)*2/kappa;
'Compare ReadIn Efficiency';
Ef_analy = 1-trapz(t,abs(Eout2).^2)/Ein_tot
Ef_analy2 = 1-(Ein_tot+4*kappa*(1/dip_amp)^2*(exp(-1/kappa*trapz(t,dipole(t).^2))-1)-2*kappa*(1/dip_amp)^2*(exp(-2/kappa*trapz(t,dipole(t).^2))-1))/Ein_tot

Ef_analy_approx = 2*kappa*(1/dip_amp)^2*(1+exp(-2/kappa*trapz(t,dipole(t).^2))-2*exp(-1/kappa*trapz(t,dipole(t).^2)))/Ein_tot

Ef_analy_approx_2 = 4*kappa*(1/dip_amp)^2*exp(-1/kappa*trapz(t,dipole(t).^2))*(cosh(1/kappa*trapz(t,dipole(t).^2))-1)/Ein_tot
test = plot_call(t,x)
%This function contains the equations to be solved.
    function dy = ydot(t,y)

	%Set the physical constants
	
	[kappa, gamma, g_0, c] = constants();
	g_0 = 1;	
	%initialize the vector to hold the solutions
        dy = zeros(2,1);
	
	%these are the equations to solve
        dy(1) = sqrt(-1)*g_0*dipole(t)*y(2)-(gamma)*y(1); %sigma

        dy(2) =  -kappa*y(2)+sqrt(2*kappa)*E_in(t) +g_0*sqrt(-1)*dipole(t)*y(1);

    end

end

%function dipole. This will be a chosen function that will be given
%for the dipole. Probably a gaussian
function out = dipole (time)
global dip_amp;
global dip_width;
[kappa, gamma, g_0 , c] = constants();
%out =-sqrt(kappa/2)/sqrt(372.2074)*exp(-(time-600).^2/(2*(210)^2))-sqrt(kappa/2)/sqrt(372.2074)*exp(-(time-2000).^2/(2*(210)^2))-sqrt(kappa/2)/sqrt(200.2074)*exp(-(time-3500).^2/(2*(210)^2));
out = dip_amp*exp(-(time-700).^2/(2*dip_width^2));
%out = dip_amp*time;
%out = exp(-(time-15).^2/20)-0.011;
%this is for width = 210, delay = 490
%out =sqrt(kappa/.005)/(17.7280)^(1/2)*exp(-(time-490).^2/(2*(10)^2));
%+kappa/2/sqrt(2*pi*(210/4)^2)*exp(-(time-2000).^2/((2*210)^2));
%out = sin( (time-120)*pi/100);
%out = time^100;
%out = time/100-time/10*sqrt(-1)

%I have calculated a funky formula for dipole that may help
%maximize our efficiency.
%N=500; %number of points to integrate over
%if length(time) ==1,
%tspan = [0:time/N:time+.001]
%out = sqrt(2/kappa*E_in(time)^2./trapz(tspan,E_in(tspan).^2)+0.001);
%else,
%out = sqrt(2/kappa*E_in(time).^2./cumtrapz(time,E_in(time).^2)+0.001);
%out = time.*0;
%end

end
%function E_in. This is the function for the incoming light
%pulse. Probably a gaussian
function out = E_in (time)
%out =  1/sqrt(2*pi*200^2)*exp(-(time-700).^2/(2*200^2));
[kappa, gamma, g_0 ,c ] =constants();
%out =   exp(-(time-700).^2/(2*200^2));
width=1000;
time;
%will make sure that a vector or 
%C. derived the optimal E_in pulse to be sent in it is given as:
%time_pulse = time(find(time<width));
%out = exp(2/kappa*time_pulse/10000)-1;
%time_pulse2 = time(find(time>=width));
%out = [out,exp(2/kappa*(width-(time_pulse2-width)))];
%time_pulse = time(find(time<width));
%will make sure that a vector or 
%C. derived the optimal E_in pulse to be sent in it is given as:
%We need to set a bound on the incoming pulse, as this will
%just keep increasing.
if length(time) > 1,
time(end)
%time_pulse = time(find(time<width));
%time_pulse2 = time(find(time>=width));
%	if time(end) > width,
	out = [exp(2/kappa*(cumtrapz(time,dipole(time).^2)))-exp(2/kappa-time*.08)]%, exp(2/kappa*(cumtrapz(time_pulse2, dipole(time_pulse2).^2)))];
%	elseif time(end) <= width,
%	out = exp(2/kappa*(cumtrapz(time_pulse,dipole(time_pulse).^2)))
%	end

length(out)
length(time)
elseif length(time) == 1 & time~=0,
time
%	if time < width,
	time = [0:time/500:time];
	out = [exp(2/kappa*(trapz(time,dipole(time).^2)))-exp(-2/kappa-time(end)*0.08)];
	
%	elseif time >= width
%	time = [width-2:time/500:time];
%	out =1-exp(2/kappa*(trapz(time, dipole(time).^2)-width));
	
%	end

elseif (length(time) == 1 & time ==0),
out = 0;
elseif length(time) < 1,
out = 0
else
'Panic'

end
%this is retrieval, so E_in = 0
out
%out = 0;
end

function [kappa, gamma, g_0, c] = constants()
kappa = 0.5;gamma = 1/200000;
c = 3.0*10^(8);
w_0 = c/765*10^(-9);
h_bar = 1.05457*10^(-34);
V = 10^(-6);
epsilon_0 = 8.854*10^(-12);
g_0 = sqrt(w_0/(2*epsilon_0*h_bar*V));
end


%This function will calculate the Delta function that represents the
%detuning present due to the external magnetic field.
function out = Delta(t)
%Phase
 %I am not sure about this, but it seems like the code is
%working on a time scale of actual nanoseconds, so I need to be consistantlength(t)

%What is N? It is the number of points to approx integrals.
%This doesn't look like it has time dependence.
%This doesn't look like it has Control. It somehow will influence
%the shape of the g_0(t) function.
Control = dipole(t)-0.00001;
C1=2.23*10^7; %gammaXe*Bx
C2=1.89*10^7; %gammaXg*BX
C3=7.52*10^7; %gammaYe ( in tesla)
C4=55.96*10^7;% gammaYg ( in tesla)
%It seems I need to normalize the field. I am hoping that control is dipole
%norma = trapz((dipole(t)*10^(9)).^2);
ControlN=Control; % normalize
Control2= ControlN.^2; %normalize, max 
W0=1.6*10^7; % detuning of the laser compare
By2=sqrt(-1)*([1:1:length(t)+2]); %change to original code, aslo removed '!!!! CHECK WITH K.
  By2= ( 2*C1*C2*C3*C4 -(1-2*Control2) .^2*(C1^2*C4^2+C2^2*C3^2)+ (1-2*Control2)*(C1*C4-C2*C3) .*sqrt((1-2*Control2) .^2*(C1*C4+C2*C3)^2-4*C1*C2*C3*C4) ) ./ (2*((1-2*Control2) .^2-1)*C3^2*C4^2); % By in function of OmegA
Delta=(0.5*sqrt(C2^2+C4^2*By2)-0.5*sqrt(C1^2+C3^2*By2)); %detuning
if numel(Delta) ~=0,
Delta=Delta - Delta(1)-4.33e+6;% laser center a the original spliting
else
Delta = 0;
end
%Recursion above? Check with K.
%Xi=D\Delta;
out = Delta*10^(-9); %This will bring Delta back into our regime (we
%work explicitlly in nanoseconds and gigahertz.
end

function out = sig_analy2(time);
global initial_cond;
[kappa, gamma, g_0, c] = constants();
steps = 500; %Number of steps in integration

%calculate the first term within the inner integral
Y1 = [0:time/steps:time; (gamma+(dipole(0:time/steps:time).^2)/kappa)];

%This will create a vector that contains the cumulative
%integration of the values of Y1 at each step
Y2 = cumtrapz(Y1(1,:), Y1(2,:));

%Now, we can caluclate the term that is the under the primary
%integral.
Y3 = sqrt(-2/kappa).*dipole(0:time/steps:time).*E_in(0:time/steps:time).*exp(Y2);
RHS =trapz(Y1(1,:),Y3);%integrate it with respect to time (Y1)

%Now, calculate the lefthand side
Y4 = [0:time/steps:time;-(gamma+(dipole(0:time/steps:time).^2)/kappa)];
out=exp(trapz(Y4(1,:), Y4(2,:)))*RHS+initial_cond(1)*exp(trapz(Y4(1,:),Y4(2,:)));
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

%I will now be constructing functions that will be examining the
%readout effifiency.
function out = analy_efout(time);
[kappa, gamma, g_0, c] = constants();
exp_arg= cumtrapz(time,-2/kappa*dipole(time).^2);
out = trapz(time, 2/kappa*dipole(time).^2.*exp(exp_arg).*exp(-2*gamma*time))
end

function out = analy_ef2out(time);
[kappa,gamma,g_0,c] = constants();
out = 1-exp(-2/kappa*trapz(time,dipole(time).^2));
end
function out = plot_call(t,x)
[kappa, gamma, g_0, c] = constants();
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
  Eout = -E_in( t ) + sqrt(2*kappa)*x(:,2)';
  Einplt = E_in( t );
  Dipolplt = dipole(t);
%for j =1:length(t);
%  analy2(j) = sig_analy2(t(j));
%  Eout_an(j) = Elec_out_analy(t(j));
%  Elec_an(j) = Elec_analy(t(j));
%  Elec_an2(j)= Elec_analy2(t(j));
%end

%This plot will investigate comparing a theoretical Eout
figure(19)
cla
plot(t,Eout,t,Einplt.*(-1+2*exp(cumtrapz(t,-Dipolplt.^2/kappa))))
h190= legend('Numeric Eout','Analytic Eout')
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
plot(t,Einplt,'--' ,t, Dipolplt,'-.', t, Eout )
h2 = legend('E_{in}', 'N^{1/2}g_0(t) GHz', 'E_{out}',1)
%set(h2,'fontsize',13);
xlabel('ns', 'fontsize',14)
ylabel('Amplitude', 'fontsize',14)
xlim([0,2000])
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

out = 1
end
