function[out] = parameters(name, num)

%parameters will contain all the modifable functions and constants 
%that can be used in the simulation of the 2-level cd memory in a cavity

	if strcmp(name,'E_in')==1,
		out = E_in(num);
	elseif strcmp(name,'constants')==1,
		out = constants(1);
	end

end

%E_in

function[out1]= E_in(input)
	x = input;
	
	out1 = x.*exp(-(x-1000)/300); %This will be normalized in the main code
end

function[out2] = constants(input)
	kappa = .5;gamma = 1/2000000000000000000000000000000000000000000;
	tau =3.4; %This will set the effective interaction time of system 
	trans_time = 3000; %this is the time when the process of reading begins
	out2 =[kappa, gamma, tau, trans_time];
end

