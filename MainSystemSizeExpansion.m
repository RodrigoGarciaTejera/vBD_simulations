%%-------------------------------------------------------------------------
% Rodrigo García-Tejera
% Contact info: rodrigo.garcia@ed.ac.uk ; rgarcia@fisica.edu.uy
% Affiliations: 
% Centre for Regenerative Medicine, University of Edinburgh, UK. 
% Facultad de Ciencias, UdelaR, Uruguay.
% Created: 19-10-2021
% Updated: 19-10-2021
%%-------------------------------------------------------------------------

%% SYSTEM ARCHITECTURE: vBD IN DIMENSIONLESS UNITS

carryingCapacity=100;       
meanConcentration=0.03;

%Propensity vectors
birthPropensity=[((1:carryingCapacity) -1).* ...
(carryingCapacity-(1:carryingCapacity)+1)/carryingCapacity /(1-meanConcentration),0];                                        
deathPropensity=[((1:carryingCapacity)),0];                                                             

%Calls constructor for BirthDeathProcess object
vBD=BirthDeathProcess(carryingCapacity,birthPropensity,deathPropensity); 

%% INITIAL CONDITIONS

%To reduce transient time initial distribution is set as a delta around the
%rate equation's nontrivial attractor 

initialDistribution=zeros(1,carryingCapacity+1);
initialDistribution(round(carryingCapacity*meanConcentration)+1)=1;

%% NUMERICALLY CALCULATED SPATIAL PROBABILITY
                                                   
[globalExtinctionState,extinctionTime]=vBD.globalExtinctionState;

%Spatial probability in half of the extinction time ensures that transient
%is over
targetDistribution=vBD.getSpatialProfile(0.5*extinctionTime,initialDistribution);

%Eliminates the extinction state and normalizes the resulting distribution
targetDistribution=targetDistribution(2:end)/sum(targetDistribution(2:end));
momentsTargetDistribution=vBD.getSpatialMoments(3,targetDistribution);


%% SYSTEM SIZE EXPANSION

%Calculates the high-order derivatives of the rate equations for the given 
%mean concentration. This function is specific for the vBD process.
dArray=dArrayValuesvBD(meanConcentration);

%Calls constructos for SSE object
expandedSystem = SystemSizeExpansion(meanConcentration,dArray); 

%Order 1 expansion
orderOneDistribution = expandedSystem.approximateDistribution(carryingCapacity,1);
momentsOrderOneDistribution=vBD.getSpatialMoments(3,orderOneDistribution(2:end)');

%Order 2 expansion
orderTwoDistribution = expandedSystem.approximateDistribution(carryingCapacity,2);
momentsOrderTwoDistribution=vBD.getSpatialMoments(3,orderTwoDistribution(2:end)');

%Order 3 expansion
orderThreeDistribution = expandedSystem.approximateDistribution(carryingCapacity,3);
momentsOrderThreeDistribution=vBD.getSpatialMoments(3,orderThreeDistribution(2:end)');


%% PLOTTING

%This can be commented after one run
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex'); 


close all
figure, pretty.plot(1:carryingCapacity,targetDistribution,'k','LineWidth',1.5), hold on
plot(0:carryingCapacity,orderOneDistribution,'r','LineWidth',1.5)
plot(0:carryingCapacity,orderTwoDistribution,'g','LineWidth',1.5)
plot(0:carryingCapacity,orderThreeDistribution,'b','LineWidth',1.5)
pretty.xylabel('Occupation number, n','Probability')
legend('Exact','$N^0$','$N^{-1/2}$','$N^{-1}$')
title(['$\phi^* =\;$',num2str(meanConcentration)],'Interpreter','Latex')

%savefig(['example_phi_0',num2str(100*meanConcentration),'_',date])
%export_fig example -pdf -transparent

% figure, pretty.plot([1:3],sqrt([(momentsOrderOneDistribution(1)-momentsTargetDistribution(1))^2, ... 
% (momentsOrderTwoDistribution(1)-momentsTargetDistribution(1))^2,...
% (momentsOrderThreeDistribution(1)-momentsTargetDistribution(1))^2]),'k','LineWidth',1.5);
% 
% figure, pretty.plot([1:3],sqrt([(momentsOrderOneDistribution(2)-momentsTargetDistribution(2))^2, ... 
% (momentsOrderTwoDistribution(2)-momentsTargetDistribution(2))^2,...
% (momentsOrderThreeDistribution(2)-momentsTargetDistribution(2))^2]),'r','LineWidth',1.5);
% 
% figure, pretty.plot([1:3],sqrt([(momentsOrderOneDistribution(3)-momentsTargetDistribution(3))^2, ... 
% (momentsOrderTwoDistribution(3)-momentsTargetDistribution(3))^2,...
% (momentsOrderThreeDistribution(3)-momentsTargetDistribution(3))^2]),'g','LineWidth',1.5);