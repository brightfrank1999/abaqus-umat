%%% Author: Fangda Cui
%%% Date: May, 2023 
%%% Neo-Hookean Material Verification
clc;
clear;


%%--Define global parameters--%%
Mu = 1.0e6;  %Shear modulus of initial phase


%%--True strain vs. True stress--%%
figure(1)
hold on;
grid on;


%%--Loading--%%
lamda = 1:0.01:3;
T1 = Mu*(lamda.^2 - 1./lamda);
plot((lamda-1)*100,T1/Mu,'k');


%%--Simulation Results--%%
[strain stress] = textread('EngineeringStrainVsTrueStress.txt','%f%f','headerlines',3);
plot(strain*100,stress/Mu,'bo');


axis([0 240 0 12]);
box on;
xlabel('Engineering Strain (%)');
ylabel('Normalized Cauchy Stress');
set(gca,'FontSize',12);
legend('Analytical Results','Simulation Results');
legend boxoff;
hold off;
