function [ t,x ,Ef, Ef2,tau,tau_list] = call_plot_new(amp,width)

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
[kappa, gamma, tau_p, trans_time,start_time] = constants();
global initial_cond;
initial_cond=y0;
global dip_amp;
dip_amp = amp;
global dip_width;
dip_width = width;
global E_norm;
global E_max;
E_max = 0;
E_norm =1;
E_norm =sqrt(trapz(E_in(0:5000).^2)); %normalize the incoming field
E_max = max(E_in(0:5000).^2); %for use in a kluge for effective fields-see plot_call
%E_norm=1;
%I want to test to see if there are any differences between the numerical normalization constant and the one achieved by analytical means.
trapz(E_in(0:5000).^2)
%Call the solver- ydot has the system of equations.
[t_out,x_out] = ode45(@ydot, tspan, y0);
 t = t_out;
x = x_out;



trapz(1/kappa*E_in(t).^2)

%To be consistant with the rest of the program, take transpose
%to get row vectors of the above quantities
Tot_Eout = 0;
Tot_Ein = 0;
%Now we want to calculate the readout efficiency.
%Hopefully, as we normalized the sigma, this should just be the
%the absolute integral of Eout squared.
Eout = sqrt(2*kappa)*x(:,2)-E_in(t);

E_in_tot_sqrd = trapz(t, E_in(t).^2);
indx_out = find(t>trans_time);
E_readout_out = sqrt(2*kappa)*x(indx_out,2) -E_in(t(indx_out));
%figure(100)
%cla
%plot(t,Eout, t, Eout2)


%Redfine the efficiency to be the maximum value of the normailzed
%population.
%Ef = max( (x(:,1).*conj(x(:,1))/E_in_tot_sqrd) );

Ef = max( (x(:,1).*conj(x(:,1))) );
Ef2 = 1-trapz(t, Eout.^2)/E_in_tot_sqrd

%Ef2 = trapz(t, Eout.^2);
tau = trapz(t(find(t<trans_time)), dipole(t(find(t<trans_time))).^2/kappa)

tau_w = trapz(t(find(t>=trans_time)), dipole(t(find(t>=trans_time))).^2/kappa)
tau_list = plot_call(t,x);
    function dy = ydot(t,y)

	%Set the physical constants
	
	[kappa, gamma, tau, trans_time, start_time] = constants();
	g_0 = 1;	
	%initialize the vector to hold the solutions
        dy = zeros(2,1);
	
	%these are the equations to solve
        dy(1) = sqrt(-1)*g_0*dipole(t)*y(2)-(gamma)*y(1); %sigma

        dy(2) =  -kappa*y(2)+sqrt(2*kappa)*E_in(t) +g_0*sqrt(-1)*dipole(t)*y(1);

    end;

end

%function dipole. This will be a chosen function that will be given
%for the dipole. Probably a gaussian
function out = dipole (time)
global dip_amp;
global dip_width;
[kappa, gamma, tau_w, trans_time, start_time] = constants();
global E_norm;
%out =-sqrt(kappa/2)/sqrt(372.2074)*exp(-(time-600).^2/(2*(210)^2))-sqrt(kappa/2)/sqrt(372.2074)*exp(-(time-2000).^2/(2*(210)^2))-sqrt(kappa/2)/sqrt(200.2074)*exp(-(time-3500).^2/(2*(210)^2));
%out = dip_amp*exp(-(time-600).^2/(2*dip_width^2));
%out = dip_amp*time;
%out = exp(-(time-15).^2/20)-0.011;
%this is for width = 210, delay = 490
%out =sqrt(kappa/.005)/(17.7280)^(1/2)*exp(-(time-490).^2/(2*(10)^2));
%+kappa/2/sqrt(2*pi*(210/4)^2)*exp(-(time-2000).^2/((2*210)^2));
%out = sin( (time-120)*pi/100);
%out = time^100;
%out = time/100-time/10*sqrt(-1)
if length(time) == 1,
if (time == 0 ) ,
	out = 0;
	return
end

	if time < trans_time; %trans time signifies end of in, begin of out
		t = (0:time/900:time);
		out =sqrt( kappa/2*(exp(2*tau_w)-1)*E_in(t(end)).^2./( (exp(2*tau_w)-1)*trapz(t,E_in(t).^2)+1));

	else %this is now readout
		t = (trans_time:(time-trans_time)/900:time);
		out =sqrt(kappa/2*E_in(t(end)-trans_time).^2./(trapz(t,-E_in(t-trans_time).^2)+1));
	%out = 0;
	end
%length(time)
else
	t_read = time(find(time<trans_time));
	t_write = time(find(time>=trans_time));
	out_read = (kappa/2 *(exp(2*tau_w)-1)*E_in(t_read).^2./(cumtrapz(t_read,(exp(2*tau_w)-1)*E_in(t_read).^2)+1)).^(.5);
	if (length(t_write)==1),
	out_write = 0;
	paste = 1
	else,
	paste = 0;
	out_write = (kappa/2*E_in(t_write-trans_time).^2./(cumtrapz(t_write,-E_in(t_write-trans_time).^2)+1)).^(.5);
	end%kluge, for problems of cumtrapz(0,0)
	out = [out_read;out_write];
%The whole purpose of trans_time is to allow me to write into the same mode that was
%written in. I want to access the original E_in function, but to do that I need to fool the
%current time into thinking that it is happening earlier. Not the most elegant way to do this, I will need to keep thinking about a better way to set up this code.

end
end
function out = dipole2(time);
global dip_amp;
global dip_width;
[kappa, gamma, tau, trans_time, start_time] = constants();
%
time_out = time;
time_out(find(time_out<trans_time))=trans_time;%heaviside on readout
%out = -dip_amp*exp(-(time-3000).^2/(2*dip_width^2))+dip_amp*exp(-(time-590).^2/(2*dip_width^2));
time_out =fliplr(time);
%out =dip_amp*exp(-(time-590).^2/(2*dip_width^2));
out = dip_amp*(exp(time/140)-1).*exp(-time/90) + dip_amp*(exp( (time_out)/140 ) -1).*exp( -(time_out)/90 );
end
function out = E_in(time)
%function E_in. This is the function for the incoming light
%out =  1/sqrt(2*pi*200^2)*exp(-(time-700).^2/(2*200^2));
[kappa, gamma, tau, trans_time, start_time ] =constants();
global E_norm;
global E_max
%This is a symmetric function that I started with
%out =   exp(-(time-900).^2/(2*300^2))/E_norm;
%time(find(time<0))=0; %this ensures that when dipole calls E_in, with the kluge of -3000, that you won't get a double read as

%Now, I will try an asymmetric function to try to illumiate the relationship
%of effective time and the symmetries

out = (time).^3.*exp(-(time-900)/400).*exp(-(time-400-900)/100)/E_norm;
indx_E_0=find(out < E_max*10^(-7)); 
out(indx_E_0) = 0.*indx_E_0;

	

end




function out =E_in2(time)
[kappa, gamma, g_0,c, start_time] = constants();
if length(time) == 1,
        if time == 0,
        out = 0;
        return  
      end

if time<2000,
time = (0:time/100:time);
out = dipole2(time(end))*exp(1/kappa*trapz(time, dipole2(time).^2));
else
out = 0; %kluge for read-out
end

else
%length(time);

out = dipole2(time).*exp(1/kappa*(cumtrapz(time,dipole2(time).^2)));
out(find(time>2000)) = 0;
end
end
function [kappa, gamma, tau, trans_time, start_time] = constants()
kappa = .5;gamma = 1/2000000000000000000000000000000000000000000;
tau = 3.4;
trans_time = 2000;
start_time = 400;
end


function out = sig_analy2(time);
global initial_cond;
[kappa, gamma, tau, trans_time, start_time] = constants();
steps = 500; %Number of steps in integration
if time == 0,
out = 0;
return
end
time_list = (0:time/steps:time)';
%Take the transpose, so that dim time_list == dim t, where
%t is the solved for time (column vector from ode45

%calculate the first term within the inner integral
Y1 = [time_list, (gamma+(dipole(time_list).^2)/kappa)];

%This will create a vector that contains the cumulative
%integration of the values of Y1 at each step
Y2 = cumtrapz(Y1(:,1), Y1(:,2));

%Now, we can caluclate the term that is the under the primary
%integral.
Y3 = sqrt(-2/kappa).*dipole(time_list).*E_in(time_list).*exp(Y2);
RHS =trapz(Y1(:,1),Y3);%integrate it with respect to time (Y1)

%Now, calculate the lefthand side
Y4 = [time_list,-(gamma+(dipole(time_list).^2)/kappa)];
out=exp(trapz(Y4(:,1), Y4(:,2)))*RHS+initial_cond(1)*exp(trapz(Y4(:,1),Y4(:,2)));
end

%this function is within the approximation of adiabatic elimantion.
function out = Elec_analy(time);
[kappa, gamma, g_0, c, start_time] = constants();
out = sqrt(-1)/kappa*dipole(time).*sig_analy2(time)+sqrt(2/kappa)*E_in(time);
end

function out =  Elec_out_analy(time);
[kappa, gamma, g_0, c, start_time] = constants();
out = sqrt(2*kappa)*Elec_analy(time)-E_in(time);
end

%this function solves the E without the adiabatic approx
function out = Elec_analy2(time,sigma);
%I don't think we can vecorize, as we are calling sig_analy2, which is
%vecorized, and passing a vector to a vecoriztion won't work well I 
%don't think. Will use loops- be slow, but should work.

[kappa, gamma, g_0, c, start_time] = constants();
%number of steps in integration
steps = 800; 
time_list = (0:time/steps:time)'; %make consistant with ode t
if time == 0,
out = 0;
return
end
Y1 = [time_list, exp(kappa*(time_list)).*(sqrt(2*kappa)*E_in(time_list)+sqrt(-1)*dipole(time_list).*sigma)];

sig_analy2(time_list)
out =1/2*exp(-kappa*time)*trapz(Y1(:,1), Y1(:,2));

end

%I will now be constructing functions that will be examining the
%readout effifiency.
function out = analy_efout(time);
[kappa, gamma, g_0, c, start_time] = constants();
exp_arg= cumtrapz(time,-2/kappa*dipole(time).^2);
out = trapz(time, 2/kappa*dipole(time).^2.*exp(exp_arg).*exp(-2*gamma*time))
end

function out = analy_ef2out(time);
[kappa,gamma,g_0,c,start_time] = constants();
out = 1-exp(-2/kappa*trapz(time,dipole2(time).^2));
end
function out = plot_call(t,x)
[kappa, gamma, tau_p, trans_time,start_time] = constants();
%Start generating the quantities to be plotted.
Eout = zeros( [length(t), 1]);
Einplt = zeros( [length(t),1] );
Dipolplt = zeros ([length(t),1] );
analy = zeros ([length(t),1]);
analy2 = zeros([length(t), 1]);
%now create plotting data sets
Eout_an = zeros( [length(t), 1]);
Elec_an = zeros( [length(t), 1]);
Elec_an2 = zeros( [length(t), 1]);
  Eout = -E_in( t ) + sqrt(2*kappa)*x(:,2);
  Einplt = E_in( t );
  Dipolplt = dipole(t);
Ein_tot = trapz(t, E_in(t).^2);
%analy2 = sig_analy2(t);
for j =1:length(t);
  analy2(j) = sig_analy2(t(j));
  Eout_an(j) = Elec_out_analy(t(j));
  Elec_an(j) = Elec_analy(t(j));
 % Elec_an2(j)= Elec_analy2(t(j));
end

%  Eout_analy = -E_in( t ) + sqrt(2*kappa)*Elec_an2(:);
%Elec_an = Elec_analy2(t);

%This plot will investigate comparing a th '.',eoretical Eout
%plot(t,analy(j))
%Now, create the plots of interest.
%sigma will be stored within x(:,1)
%E will be stored within x(:,2) 
%the (:,1) notation means 'all the rows in colulmn 1'
tplt = t;
figure(1)
clf
%set(1,'defaultaxesfontsize',15)
plot(t,x(:,1),t,imag(x(:,1)) ,t, x(:,2),t,imag(x(:,2)))
 legend('Re \sigma','Im \sigma','Re E','Im E',1)
xlabel('ns','fontsize',14)
ylabel('Amplitude','fontsize',14)
%set(h,'fontsize',13)
%print('-dpdf', strcat('sige',num2str(dip_amp),'.pdf'))


figure(3)
clf
plot(t,x(:,1).*conj(x(:,1)),'--', t,x(:,2).*conj(x(:,2)));
 legend('|\sigma|^2', '|E|^2')
xlabel('ns','fontsize',14)
ylabel('Amplitude','fontsize',14)
%print('-dpdf',strcat('badrunexcite',num2str(dip_amp),'.pdf'))


figure(2)
cla
clf
plot(t,abs(Einplt),'--' ,t,abs( Dipolplt),'-.', t, abs(Eout) )
h2 = legend('E_{in}', 'N^{1/2}g_0(t) GHz', 'E_{out}',1)
%set(h2,'fontsize',13);
xlabel('ns', 'fontsize',14)
ylabel('Amplitude', 'fontsize',14)
%xlim([0,2000])
%
figure(101)
clf 
cla
plot(t, x(:,1).*conj(x(:,1))/Ein_tot,t, x(:,2).*conj(x(:,2))/Ein_tot);
legend('|\sigma|^2 / |E_{in}|^2','|E|^2 / |E_{in}|^2')
xlabel('ns')
ylabel('Normalized Population')




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
indx_read =find(t>=trans_time);
indx_write = find(t<trans_time);
effective_write = sqrt(kappa)*Einplt./(Dipolplt+eps);
Eout_max = max(Eout);

global E_max;
indx_Eout=find(Dipolplt <E_max*10^(-4));
effective_read=-sqrt(kappa)*Eout./(Dipolplt+eps);
yu = effective_read;
effective_read(indx_Eout)=0*indx_Eout;
indx_Eout
yu = Eout
yu(indx_Eout) = 666
[Eout,Dipolplt+eps,yu]
effective_read(indx_write)=0;
effective_write(indx_read) =0;

%Now, we need to smooth the data to get rid of /0 noise
span = 4;
window = ones(span,1)/span;
effective_write = conv(effective_write, window, 'same');
effective_read = conv(effective_read, window, 'same');
global E_norm;
figure(123)
cla
clf
effective_fields = effective_write+effective_read;
subplot(211)
plot(t, effective_write, t, effective_read)
h12 = legend('Effective Write', 'Effective Read',1)
title('Effective Fields')
xlabel('Time')
ylabel('Intensity')
subplot(212)
E_norm
ef_analy_write= sqrt(2/(exp(2*3.4)-1))*exp(cumtrapz(t, Dipolplt.^2)/kappa);
ef_analy_read= [indx_write*0;sqrt(2)*exp(-cumtrapz(t(indx_read), Dipolplt(indx_read).^2)/kappa)];
ef_analy_write(indx_read)=0;
ef_analy_read(indx_write)=0;
plot(t, ef_analy_write,t , ef_analy_read);
h13= legend('Effective Write', 'Effective Read','1')
figure(190)
clf
subplot(311)
plot(t,abs(Einplt),'--' ,t, abs(Eout), '.',t, abs(Eout_an) )
h2 = legend('E_{in}', 'E_{out}', 'Analytic E_{out}',1)
%set(h2,'fontsize',13);
xlabel('ns', 'fontsize',14)
ylabel('Amplitude', 'fontsize',14)
xlim([0,3500])
%
subplot(312)
plot(tplt,abs( Dipolplt),'-.')
h3 = legend('g(t)')
xlabel('ns')
ylabel('Frequency (GHz)')
xlim([0,3500])
subplot(313) 
tau_list = cumtrapz(t,Dipolplt.^2/kappa);
tau_list_w = tau_list(find(t<trans_time));
tau_list_r = tau_list(find(t>= trans_time));
chx2=plot(cumtrapz(t,Dipolplt.^2/kappa), effective_write, cumtrapz(t,Dipolplt.^2/kappa), effective_read, tau_list_w, sqrt(2)*exp(tau_list_w)/max(exp(tau_list_w)),tau_list_r, sqrt(2)*exp(-tau_list_r)/max(exp(-tau_list_r)))
h12 = legend('Effective Write', 'Effective Read',1)
xlabel('Effective Time')
ylabel('Intensity')
xlim([0,2*3.4])

out = tau_list;
end
