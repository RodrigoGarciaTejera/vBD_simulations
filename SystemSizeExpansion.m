%%-------------------------------------------------------------------------
% Rodrigo García-Tejera
% Contact info: rodrigo.garcia@ed.ac.uk ; rgarcia@fisica.edu.uy
% Affiliations: 
% Centre for Regenerative Medicine, University of Edinburgh, UK. 
% Facultad de Ciencias, UdelaR, Uruguay.
% Created: 19-10-2021
% Updated: 19-10-2021
%%-------------------------------------------------------------------------

%This class generates and handles the expanded approximate solutions of a 
%master euqation, obtained from the system size expansion. Each object 
%represents a master equation, its attributes being the Jacobian of the 
%rate equations, the dispersion obtained by the linear noise approximation, 
%and a matrix D composed of combinations of the high-order derivatives of
%the propensity functions. This class can only deal with ME from systems 
%that have a single asymptotically stable fixed point, and propensities 
%that are proportional to the system size.   

% Order of expansion coefficients:order 2: a11-a31 order 3: a22-a42-a62%
% Order for dArray:order 1:D11(Jacobian)D20 (related to variance) 
% order 2: D12-D21-D30 order 3:D22-D13-D40-D31

classdef SystemSizeExpansion 
    properties
        dArray
        jacobian
        variance
        meanConcentration
    end
    
    methods
        
        %Constructor for the class's objects
        function obj = SystemSizeExpansion(meanConcentration,dArray)            
            obj.dArray=dArray;
            obj.jacobian=dArray{1};
            obj.meanConcentration=meanConcentration;
            obj.variance=-dArray{2}/2/obj.jacobian;
        end     
        
        %Returns the coefficients of the approximate solution inspired in
        %the SSE up to a desired order, for the SSE instance given. Order 
        %must be higher than one.
        function expansionCoefficients = getCoefficients(obj,order)
              
            if order==2 
                expansionCoefficients=cell(1,2);
                
                expansionCoefficients{1}=-obj.variance*obj.dArray{3}/2/obj.jacobian;
                
                expansionCoefficients{2}=-obj.variance^2*obj.dArray{3}/6/obj.jacobian - ... 
                obj.variance*obj.dArray{4}/6/obj.jacobian -obj.dArray{5}/18/obj.jacobian;            
            
            elseif order==3
                expansionCoefficients=obj.getCoefficients(order-1);
                
                expansionCoefficients{3}=-expansionCoefficients{1}*(obj.dArray{4}/4/obj.jacobian + ... 
                3 * obj.variance *obj.dArray{3}/4/obj.jacobian)- ... 
                expansionCoefficients{2}*3*obj.dArray{3}/2/obj.jacobian-obj.variance*obj.dArray{6}/8/obj.jacobian - ... 
                obj.variance^2*obj.dArray{7}/4/obj.jacobian;
            
                expansionCoefficients{4}=-expansionCoefficients{1}*(obj.dArray{5}/24/obj.jacobian + ...
                obj.variance*obj.dArray{4}/8/obj.jacobian + obj.variance^2*obj.dArray{3}/8/obj.jacobian)- ...
                obj.dArray{8}/96/obj.jacobian - obj.variance*obj.dArray{9}/24/obj.jacobian - ...
                obj.variance^2*obj.dArray{6}/16/obj.jacobian - obj.variance^3*obj.dArray{7}/24/obj.jacobian - ...
                expansionCoefficients{2}*(3*obj.dArray{4}/8/obj.jacobian + 7*obj.variance*obj.dArray{3}/8/obj.jacobian);
                
                expansionCoefficients{5}=0.5*expansionCoefficients{2}^2;
           
            elseif order==4
                expansionCoefficients=obj.getCoefficients(order-1);
                expansionCoefficients{6}=-obj.dArray{3}*expansionCoefficients{3}/obj.jacobian;
                expansionCoefficients{7}=-(expansionCoefficients{1}/3/obj.jacobian)*(3*obj.dArray{6}*obj.variance/4+obj.dArray{9}/6);
                expansionCoefficients{8}=0;                                %THESE ARE 0 IN THE vBD CASE ONLY, CORRECT THAT TO GENERALIZE 
                expansionCoefficients{9}=0;
                expansionCoefficients{10}=0;
% Order of expansion coefficients:order 2: a11-a31 order 3: a22-a42-a62%
% Order for dArray:order 1:D11(Jacobian)D20 (related to variance) 
% order 2: D12-D21-D30 order 3:D22-D13-D40-D31
            else 
                expansionCoefficients=[];
                
            end
            
        end
        
        %Approximates solution of the spatial component of the ME up to a
        %desired order of the SSE
        function proxyPDF = approximateDistribution(obj,carryingCapacity,order)
            
            meanOccupation=carryingCapacity*obj.meanConcentration;
            occupationVector=[0:carryingCapacity] - meanOccupation;
            occupationVariance=carryingCapacity*obj.variance;
            
            %Order 1 approximation (discrete version of LNA)
            proxyPDF = 0.5*exp(- occupationVector.^2 /2 /occupationVariance)/sqrt(2*pi*occupationVariance)...
            .* (erfz((1i*occupationVector+pi*occupationVariance)/(sqrt(2*occupationVariance)))- ...
            erfz((1i*occupationVector-pi*occupationVariance)/(sqrt(2*occupationVariance))));
            
            if order > 1 
                
                expansionCoefficients = obj.getCoefficients(order);
                
                firstDerivative=-occupationVector/occupationVariance.*proxyPDF+ ... 
                exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);
                
                secondDerivative=-proxyPDF/occupationVariance - occupationVector/occupationVariance .*firstDerivative + ...
                pi*exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
                
                thirdDerivative=-2/occupationVariance * firstDerivative - ... 
                occupationVector/occupationVariance .* secondDerivative - ... 
                pi^2 * exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);
                
                %Order 2 approximation
                proxyPDF= proxyPDF - expansionCoefficients{1}*firstDerivative - expansionCoefficients{2}* ...
                carryingCapacity*thirdDerivative;
                
                if order > 2
                    fourthDerivative=-3/occupationVariance * secondDerivative - occupationVector/occupationVariance .* ...
                    thirdDerivative - pi^3 * exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
                    
                    fifthDerivative=-4/occupationVariance * thirdDerivative - occupationVector/occupationVariance .* ...
                    fourthDerivative + pi^4 * exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);
                    
                    sixthDerivative=-5/occupationVariance * fourthDerivative - occupationVector/occupationVariance .* ...
                    fifthDerivative + pi^5 * exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
                     
                    %Order 3 approximation
                    proxyPDF = proxyPDF + expansionCoefficients{3}*secondDerivative+ ...
                    expansionCoefficients{4}*carryingCapacity*fourthDerivative + ...
                    expansionCoefficients{5}*carryingCapacity^2*sixthDerivative;
                
                end      
            end
            
        end


        function proxyPDF = approximateDistributionRenormalized(obj,carryingCapacity,order)
            
            
            expansionCoefficients = obj.getCoefficients(4);          
            
            switch order
                case 2 
                    occupationVariance=carryingCapacity*obj.variance;
                    meanFluctuations=expansionCoefficients{1}/sqrt(carryingCapacity);
                case 3 
                    occupationVariance=carryingCapacity*obj.variance+2*expansionCoefficients{3}-expansionCoefficients{1}^2;
                    meanFluctuations=expansionCoefficients{1}/sqrt(carryingCapacity);
                case 4
                    occupationVariance=carryingCapacity*obj.variance+2*expansionCoefficients{3}-expansionCoefficients{1}^2;
                    meanFluctuations=expansionCoefficients{1}/sqrt(carryingCapacity)+expansionCoefficients{6}/(carryingCapacity^(3/2));
            end

     

            meanOccupation=carryingCapacity*obj.meanConcentration+sqrt(carryingCapacity)*meanFluctuations;
            occupationVector=[0:carryingCapacity] - meanOccupation;
            
            %RENORMALIZATION OF COEFFIEICNTS
  
            expansionCoefficientsRen{1}=0; expansionCoefficientsRen{2}=expansionCoefficients{2};
            expansionCoefficientsRen{3}=0; 
            expansionCoefficientsRen{4}=expansionCoefficients{4}-expansionCoefficients{1}*expansionCoefficients{2};
            expansionCoefficientsRen{5}=expansionCoefficients{5};
            expansionCoefficientsRen{6}=0;
            expansionCoefficientsRen{7}=expansionCoefficients{1}*(expansionCoefficients{1}^2-2*expansionCoefficients{3})+expansionCoefficients{7};
            expansionCoefficientsRen{8}=expansionCoefficients{8}-expansionCoefficients{4}*expansionCoefficients{1};
            
            expansionCoefficientsRen{9}=expansionCoefficients{9}-expansionCoefficients{1}*expansionCoefficients{5};
            expansionCoefficientsRen{10}=expansionCoefficients{10};
           


% Order of expansion coefficients:order 2: a11-a31 order 3: a22-a42-a62%

            expansionCoefficients=expansionCoefficientsRen;

            %Order 1 approximation (discrete version of LNA)
            proxyPDF = 0.5*exp(- occupationVector.^2 /2 /occupationVariance)/sqrt(2*pi*occupationVariance)...
            .* (erfz((1i*occupationVector+pi*occupationVariance)/(sqrt(2*occupationVariance)))- ...
            erfz((1i*occupationVector-pi*occupationVariance)/(sqrt(2*occupationVariance))));

                firstDerivative=-occupationVector/occupationVariance.*proxyPDF+ ... 
                exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);
                
                secondDerivative=-proxyPDF/occupationVariance - occupationVector/occupationVariance .*firstDerivative + ...
                pi*exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
                
                thirdDerivative=-2/occupationVariance * firstDerivative - occupationVector/occupationVariance .* ...
                secondDerivative - pi^2 * exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);
                
                %Order 2 approximation
                proxyPDF= proxyPDF - expansionCoefficients{1}*firstDerivative - expansionCoefficients{2}* ...
                carryingCapacity*thirdDerivative;
              
                fourthDerivative=-3/occupationVariance * secondDerivative - occupationVector/occupationVariance .* ...
                thirdDerivative - pi^3 * exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
                    
                fifthDerivative=-4/occupationVariance * thirdDerivative - occupationVector/occupationVariance .* ...
                fourthDerivative + pi^4 * exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);
                    
                sixthDerivative=-5/occupationVariance * fourthDerivative - occupationVector/occupationVariance .* ...
                fifthDerivative + pi^5 * exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
                     
                %Order 3 approximation
                proxyPDF = proxyPDF + expansionCoefficients{3}*secondDerivative+ ...
                expansionCoefficients{4}*carryingCapacity*fourthDerivative + ...
                expansionCoefficients{5}*carryingCapacity^2*sixthDerivative;        

                seventhDerivative=-6/occupationVariance * fifthDerivative - occupationVector/occupationVariance .* ...
                sixthDerivative - pi^6 * exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);

                eighthDerivative=-7/occupationVariance * sixthDerivative - occupationVector/occupationVariance .* ...
                seventhDerivative - pi^7 * exp(-pi^2*occupationVariance/2)*cos(pi*occupationVector)/(pi*occupationVariance);
 
                ninethDerivative=-8/occupationVariance * seventhDerivative - occupationVector/occupationVariance .* ...
                eighthDerivative + pi^8 * exp(-pi^2*occupationVariance/2)*sin(pi*occupationVector)/(pi*occupationVariance);

                %Order 4 approximation
               
%                 proxyPDF = proxyPDF-carryingCapacity^(-1)*expansionCoefficients{6}*firstDerivative - ...
%                 carryingCapacity^(0)*expansionCoefficients{7}*thirdDerivative - ... 
%                 carryingCapacity^1*expansionCoefficients{8}*fifthDerivative ...
%                 - carryingCapacity^2*expansionCoefficients{9}*seventhDerivative- ... 
%                 carryingCapacity^3*expansionCoefficients{10}*ninethDerivative;

                proxyPDF=proxyPDF/sum(proxyPDF);

        end    
    end
end
