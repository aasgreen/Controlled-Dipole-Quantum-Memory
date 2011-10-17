function [ t,x ,ef_write, ef_read, tau_write, tau_read] = simple_simu_cd()

%simple_simu_cd.m will solve the coupled ODE's for sigma and E. 
%   The coupled ODE's for sigma and E are in the function dy.You can specify the
%parameters and input field in the parameters.m file that is also included.
% You can run this by entering the matlab environment and typing:
% [t,x,ef_write, ef_read, tau_write, tau_read] = simple_simu_cd() 

%Main Script to run:

%First, generate the time vector and the initial conditions: 
tspan = [0 8000];
y0 = [0,0] %retrieval, so |sigma(t=0)| = 1
[kappa, gamma, tau, trans_time] = constants();
global initial_cond;
initial_cond=y0;

%Normalize the incoming E_in first
global E_norm;
E_norm =1;
E_norm =sqrt(trapz(E_in(0:8000).^2)); %normalize the incoming field
%E_norm is a global variable, so it will be a part of the E_in function now


%Call the solver- ydot has the system of equations.
[t_out,x_out] = ode45(@ydot, tspan, y0);
 t = t_out;
x = x_out;
%Now calcualte the efficiencies of the system:
[ef_write, ef_read, tau_write, tau_read] = eff_tau(t,x)

oop = plot_call(t,x);
    function dy = ydot(t,y)

	%Set the physical constants
	
	[kappa, gamma, tau, trans_time] = constants();
	%initialize the vector to hold the solutions
        dy = zeros(2,1);
	
	%these are the equations to solve
        dy(1) = sqrt(-1)*dipole(t)*y(2)-(gamma)*y(1); %sigma

        dy(2) =  -kappa*y(2)+sqrt(2*kappa)*E_in(t) +sqrt(-1)*dipole(t)*y(1);

    end;

end

%function dipole. This will be a chosen function that will be given
%for the dipole. Probably a gaussian
function out = dipole (time)
	[kappa, gamma, tau, trans_time] = constants();
	tau_r = 3.4;
	if time == 0,
	out = 0;
	return
	end


	if time < trans_time; %trans time signifies end of in, begin of out
		t = (-eps:time/1000:time);
		out =sqrt( kappa/2*(exp(2*tau)-1)*E_in(t(end))^2/( (exp(2*tau)-1)*trapz(t,E_in(t).^2)+1));

	else %this is now readout
		if E_in(time-trans_time) < 1E-9,
			out = 0;
			return
		end
		%The above is neccesary, as the dividing by a small number
		% involved in the formula accentuates numeric error
		t = (trans_time:(time-trans_time)/100:time);
		out =sqrt(kappa/2*(1-exp(-2*tau_r))*E_in(t(end)-trans_time)^2/(trapz(t-trans_time,-(1-exp(-2*tau_r))*E_in(t-trans_time).^2)+1));
	%out = 0;
	end
end


function out = E_in(time)
[kappa, gamma, tau, trans_time ] =constants();
global E_norm;
out = parameters('E_in', time)/E_norm;

end



function [kappa, gamma, tau, trans_time] = constants()
raw = parameters('constants');
kappa =raw(1);
gamma = raw(2);
tau = raw(3);
trans_time = raw(4);
end

function [ef_write, ef_read, tau_write, tau_read] = eff_tau(t,x)
%eff_tau will calculate the efficiens of the various processes, as 
%well as displaying the tau (effective time).
[kappa, gamma, tau, trans_time] = constants()
E_in_list = E_in(t);
Eout = sqrt(2*kappa)*x(:,2)-E_in_list;
dipole_list = zeros([length(t),1]);
for j=1:length(t);
dipole_list(j) = dipole(t(j));
end
E_in_tot_sqrd = trapz(t, E_in(t).^2);
indx_read = find(t>=trans_time);
indx_write = find(t<trans_time);

%The write efficiency is calculated to be the maximum value of the
%normalized atomic population. Some care needs to be exersized with this, as 
%sometimes, during the process of writing the photon to the atomic ensemble, there
%is a rapid spike, or local maximum that can throw this measurment off. However,
%an examination of the plot of the atomic population is a good check of this.
ef_write = max( (x(:,1).*conj(x(:,1))) );

%The read efficiency is calculated as the integral of the read field, after
%the trans_time. This ensures that no 'leaky field' that occured during
%the write process will influence the read efficiency. Ensure that the incoming
%field is normalized, or this calculation won't be accurate. Also, under certain
%circumstances the code can become pathological and the E_out can tend
%to infinity, usually when large discontinunities are present. Just something to
%watch for.
ef_read = trapz(t(indx_read),Eout(indx_read).^2);
indx_write
tau_write = trapz(t(indx_write), dipole_list(indx_write).^2/kappa);

tau_read = trapz(t(indx_read), dipole_list(indx_read).^2/kappa);
end


function out = plot_call(t,x)
%This function will create all the plots of interest.
[kappa, gamma, tau, trans_time] =constants();

%Start generating the quantities to be plotted.
Eout = zeros( [length(t), 1]);
Einplt = zeros( [length(t),1] );
Dipolplt = zeros ([length(t),1] );

Einplt = E_in(t)
Eout = -Einplt + sqrt(2*kappa)*x(:,2);

%To get rid of some of the complexity, I removed the vectorization of 
%the dipole function, a for loop is now neccesary to calculate things
for j=1:length(t);
Dipolplt(j) = dipole(t(j));
end
Ein_tot = trapz(t, Einplt.^2);


%PLOTTING

%The first plot describes the two solutions of the equations of motion
figure(1)
clf
plot(t,x(:,1),t,imag(x(:,1)) ,t, x(:,2),t,imag(x(:,2)))
 legend('Re \sigma','Im \sigma','Re E','Im E',1)
xlabel('ns','fontsize',14)
ylabel('Amplitude','fontsize',14)
%set(h,'fontsize',13)
%print('-dpdf', strcat('sige',num2str(dip_amp),'.pdf'))

%PLOT 2: Describes the square of these terms
figure(2)
clf
plot(t,x(:,1).*conj(x(:,1)),'--', t,x(:,2).*conj(x(:,2)));
 legend('|\sigma|^2', '|E|^2')
xlabel('ns','fontsize',14)
ylabel('Amplitude','fontsize',14)
%print('-dpdf',strcat('badrunexcite',num2str(dip_amp),'.pdf'))

%PLOT 3: Describes the interaction between the fields and the coupling
figure(3)
cla
clf
plot(t,abs(Einplt),'--' ,t,abs( Dipolplt),'-.', t, abs(Eout) )
h2 = legend('E_{in}', 'N^{1/2}g_0(t) GHz', 'E_{out}',1)
%set(h2,'fontsize',13);
xlabel('ns', 'fontsize',14)
ylabel('Amplitude', 'fontsize',14)
%xlim([0,2000])
%

%PLOT 4: If the E_in field is not normalized, this plot will give the normalized
%results of the abs square of the solutions
figure(4)
clf 
cla
plot(t, x(:,1).*conj(x(:,1))/Ein_tot,t, x(:,2).*conj(x(:,2))/Ein_tot);
legend('|\sigma|^2 / |E_{in}|^2','|E|^2 / |E_{in}|^2')
xlabel('ns')
ylabel('Normalized Population')



%This next series of plots is concerned with comparing various effective
%quantities with their 'physical' components.
indx_read =find(t>=trans_time);
indx_write = find(t<trans_time);
effective_write = sqrt(kappa)*Einplt./(Dipolplt+eps);
effective_read=-sqrt(kappa)*Eout./(Dipolplt+eps);

%This next step is to 'smooth' out some values. If it is not included, then
%a lot of noise crops up in times sooner that trans_time. This is probably due to the fact that above, the effective fields involve dividing two very small numbers.
effective_read(indx_write)=0;


%PLOT 5: This plot contains two plots. The first shows the numerically calculated
%effective fields, based on what the dipole and the Ein or Eout fields were. The 
%second plot shows the solution of the equations as a function of time. This would%be the analytical solution that I am comparing to.
figure(5)
cla
clf
effective_fields = effective_write+effective_read;
subplot(211)
plot(t, effective_write, t, effective_read)
h12 = legend('Effective Write (Numeric)', 'Effective Read (Numeric)',1)
title('Effective Fields')
xlabel('Time')
ylabel('Amplitude')
subplot(212)
ef_analy_write= sqrt(2/(exp(2*tau)-1))*exp(cumtrapz(t, Dipolplt.^2)/kappa);
ef_analy_read= [indx_write*0;sqrt(2)*exp(-cumtrapz(t(indx_read), Dipolplt(indx_read).^2)/kappa)];
ef_analy_write(indx_read)=0;
ef_analy_read(indx_write)=0;
plot(t, ef_analy_write, t, ef_analy_read);
h13= legend('Effective Write (Analytic)', 'Effective Read (Analytic)','1')
xlabel('Time')
ylabel('Amplitude')



%one last thing to try:
effec_eout = [0*indx_write;sqrt(2*(1-cumtrapz(Eout(indx_read).^2)))];
%Plot 6: This plot compares the 'physical' fields to the equivalent 'effective' fields. Note, the first plot is in time, the second one is in effective time tau.
figure(6)
clf
cla
subplot(211)
plot(t,abs(Einplt),'--' ,t,abs( Dipolplt),'-.', t, abs(Eout) )
h2 = legend('E_{in}', 'N^{1/2}g_0(t) GHz', 'E_{out}',1)
%set(h2,'fontsize',13);
xlabel('ns', 'fontsize',14)
ylabel('Amplitude', 'fontsize',14)
%
subplot(212)
cla
plot(cumtrapz(t,Dipolplt.^2/kappa), effective_write, cumtrapz(t,Dipolplt.^2/kappa), effec_eout)
h12 = legend('Effective Write', 'Effective Read',1)
xlabel('\tau')
ylabel('Amplitude')


out = 1
end
