function [test_ef] = run_ef();
indx = [1: 30]
plot_width=zeros(30,30);
plot_dip =zeros(30,30);
ef_numer = zeros(30,30);
ef_analy = zeros(30,30);
test_ef = zeros(30,4);
for i=1:length(indx),	
	  dip = 0.01+(0.09-0.01)/30*i;
	
%	for j=1:length(indx),
%	 width = 140+(230-140)/30*j; 
	  [t,x,Efnum_ideal,Efnum_notideal, tau] = call_plot_new(dip,210);
%  	 plot_width(i,j) = width;
%	 plot_dip(j,i) = dip;
%	 ef_numer(i,j) = Efnum
%  	 ef_analy(i,j) = Efanaly;
%	[j,i]
%	Efnum
        	test_ef(i,1) = dip;
	test_ef(i,2) = Efnum_ideal;
	test_ef(i,3) = Efnum_notideal;
	test_ef(i,4) = tau
	i
 % 	 end
     end	

end
