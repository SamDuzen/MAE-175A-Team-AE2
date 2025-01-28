fit_var = polyfit(loc_u,Cp_u_zero_deg_20,2);
fit_x = linspace(loc(2),loc(end),1000);
fit_line = polyval(fit_var,fit_x);
figure(1)
hold on
plot(fit_x,fit_line)
plot(loc_u,Cp_u_zero_deg_20)
hold off

%%
cpu020_pval = polyfit(loc_u,Cp_u_zero_deg_20,2);
fit_var = cpu020_pval;
fit_x = linspace(0,3.5,1000);
fit_line = polyval(fit_var,fit_x);
figure(1)
hold on
plot(fit_x,polyval(cpu020_pval,fit_x),'-b')
plot(loc_u,Cp_u_zero_deg_20,':b')
plot(fit_x,polyval(cpl020_pval,fit_x),'-r')
plot(loc_l,Cp_l_zero_deg_20,':r')
plot(fit_x,dum_func(fit_x),'-m')
%plot(chord_spread,test,'.g')
hold off
%legend


%%
dum_func = @(x) (cpn020_pval(1)/3*x.^3 + cpn020_pval(2)/2*x.^2 + cpn020_pval(3)*x)/3.5;

figure(2)
plot(fit_x,dum_func(fit_x))

%% Code for function creation
test = calculate_cpn_values(cpn020_pval,3.5,chord_spread);
plot(chord_spread,test)
function output_values = calculate_cpn_values(cpn,c,input_values)
calculate_values = @(x) (cpn(1)/3*x.^3 + cpn(2)/2*x.^2 + cpn(3)*x)/c;
output_values = calculate_values(input_values);
end


% dummy_variable = CP_NAME;
% dummy_function = @(x) (dummy_variable(1)/3*x.^3 + dummy_variable(2)/2*x.^2 + dummy_variable(3)*x)/3.5;
% FUNCTION_NAME = dummy_function;