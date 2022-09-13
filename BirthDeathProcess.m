%%-------------------------------------------------------------------------
% Rodrigo García-Tejera
% Contact info: rodrigo.garcia@ed.ac.uk ; rgarcia@fisica.edu.uy
% Affiliations: 
% Centre for Regenerative Medicine, University of Edinburgh, UK. 
% Facultad de Ciencias, UdelaR, Uruguay.
% Created: 18-10-2021
% Updated: 25-10-2021
%%-------------------------------------------------------------------------

%This class generates and handles one-step birth and death processes with
%carrying capacity. It does not work with unbounded birth and death 
%processes, nor with chemical reactions networks featuring changes of more
%than one molecule per reaction.


classdef BirthDeathProcess 
    
    properties
        carryingCapacity
        birthPropensity
        deathPropensity
    end

    methods
        
        %Constructor for BirthDeathProcess objects
        function obj= BirthDeathProcess(carryingCapacity,birthPropensity,deathPropensity)
            
            obj.carryingCapacity=carryingCapacity;
            obj.birthPropensity=birthPropensity;
            obj.deathPropensity=deathPropensity;
        
        end     
        
        %Returns the transition matrix given the propensities and carrying
        %capacity
        function transitionMatrix=getTransitionMatrix(obj)
      
            transitionMatrix=-diag(obj.birthPropensity+ ... 
            [0,obj.deathPropensity(1:end-1)])+diag(obj.deathPropensity(1:end-1),1)+ ...
            diag(obj.birthPropensity(1:end-1),-1);
        
        end
        
        %Returns the occupation probability distribution as a function of
        %time for the given BirthDeathProcess instance and initial 
        %distribution, calculated by numerical matrix exponentiation.       
        function probabilityDistribution=getprobabilityDistribution(obj,timeVector,initialProbability)
            
            transitionMatrix=obj.getTransitionMatrix;
            probabilityDistribution=zeros(numel(timeVector),obj.carryingCapacity+1);
            probabilityDistribution(1,:)=initialProbability;
            
            for k=2:numel(timeVector)
            probabilityDistribution(k,:)=expm(transitionMatrix*timeVector(k))*initialProbability'; 
            end
            
        end
        
        %Returns the spatial profile at a given time point for a 
        %BirthDithProcess instance, given the initial conditions.
        %obj.getSpatialProfile(obj,time,initial_probability) returns the
        %the second column of 
        %obj.getProabilityDistribution([0,time],initial_probability)
        
        function spatialProbability=getSpatialProfile(obj,time,initial_probability)
            
            transitionMatrix=obj.getTransitionMatrix;
            spatialProbability=expm(transitionMatrix*time)*initial_probability'; 
        
        end
        
        %The boolean globalExtinctionState is true if the BirthDithProcess
        %instance features a global extinction state, and if so,
        %exctinctionTime gives the expected extinction time, calculated as
        %the real part of the transition matrix's largest non-zero eigevnalue
        function [globalExtinctionState,extinctionTime]=globalExtinctionState(obj)
            
            transitionMatrix=obj.getTransitionMatrix;
            sortedEigenvalues=sort(real(eig(transitionMatrix)),'descend');
            
            toleranceForZeroEigenvalue=5e-16;
            if sortedEigenvalues(1) <= toleranceForZeroEigenvalue
                
                globalExtinctionState=true;                
                positionFirsNegativeEigenvalue=find(sortedEigenvalues<= -toleranceForZeroEigenvalue,1,'first');
                extinctionTime=-1/sortedEigenvalues(positionFirsNegativeEigenvalue);
            else
                
                globalExtinctionState=false; 
                extinctionTime=[];
            end

        end
        
        function FanoFactor=getFanoFactor(obj,timeVector, initialProbability,varargin)
            
            switch nargin 
                
                case 3
                    probabilityDistribution=obj.getprobabilityDistribution(timeVector,initialProbability);
                    
                    occupationMean=probabilityDistribution*(0:obj.carryingcapacity)';
                    occupationVariance=probabilityDistribution*((0:obj.carryingcapacity).^2)' - occupationMean.^2;
                    FanoFactor=occupationVariance./occupationMean;
                    
                case 4
                  
                excludeExtinction=varargin{:} ;
                
                if excludeExtinction == false
                    FanoFactor=obj.getFanoFactor(timeVector, initialProbability);
                
                elseif excludeExtinction == true
                    probabilityDistribution=obj.getprobabilityDistribution(timeVector,initialProbability);
                    probabilityDistribution(:,1)=[];
                    
                    
                    probabilityDistribution=probabilityDistribution./sum(probabilityDistribution,2);
                    
                    occupationMean=probabilityDistribution*(1:obj.carryingCapacity)';
                    occupationVariance=probabilityDistribution*((1:obj.carryingCapacity).^2)' - occupationMean.^2;
                    FanoFactor=occupationVariance./occupationMean;
                end
            end
        end
        
        function moments=getSpatialMoments(obj,order,spatialProbability)
            
            moments=zeros(1,order);
            occupationNumbers=1:obj.carryingCapacity;
            for r=1:order
                try
                moments(r)=sum(spatialProbability'.*occupationNumbers.^r);
                catch
                    keyboard
                end
            end   
        end

%Assumes that the system has an absorbing boundary at zero
        function tau=firstPassageTimes2Extinction(obj)
            
            tau=zeros(1,obj.carryingCapacity);                             %tau is the vector of first passage times
            tau(1)=sum(flip(cumprod([1,obj.birthPropensity(2:(end-1))])) ... 
            ./flip(cumprod(obj.deathPropensity(1:end-1))));

            for r=2:obj.carryingCapacity
            tau(r)=tau(r-1) + sum(flip(cumprod([1,obj.birthPropensity(r+1:(end-1))])) ... 
            ./flip(cumprod(obj.deathPropensity(r:(end-1)))));
            end
       
       %if isnan(sum(tau))
       %    keyboard
       %end
             
        end

%Assumes that the system has an absorbing boundary at zero
% Modes: 1-->Making the rate of the transition 1-->0 and P0 be zero.
%        2-->Balancing P1 with the inclusion of the transition rate 1-->0
%            but making P0=0. 

        function QSS=getGardinerQss(obj,mode)

            if mode==1

                factor1=prod(obj.deathPropensity(2:end-1)./obj.birthPropensity(2:end-1));
                factor2=sum(cumprod(flip(obj.deathPropensity(3:end-1))./flip(obj.birthPropensity(3:end-1))));
                normalizationFactor=1/(1+factor1+factor2);
                
            elseif mode==2

                factor1=prod(obj.deathPropensity(2:end-1)./[obj.birthPropensity(2)+ ... 
                obj.deathPropensity(1),obj.birthPropensity(3:end-1)]);
                factor2=sum(cumprod(flip(obj.deathPropensity(3:end-1))./flip(obj.birthPropensity(3:end-1))));
                normalizationFactor=1/(1+factor1+factor2);
            end

            QSS=[0,factor1*normalizationFactor, flip(cumprod(flip(obj.deathPropensity(3:end-1)) ... 
            ./flip(obj.birthPropensity(3:end-1))))* normalizationFactor,normalizationFactor];
            
        end

        function PWKB=getvBDWKB(obj,meanConcentration)
        N=obj.carryingCapacity;
        phi=meanConcentration;
        x=(1:(N-1))/N;

        PWKB=sqrt(N./(2*pi*(1-x))).*(phi./x).*exp(N*(phi-x)).*((1-x)./(1-phi)).^(N*(x-1));
        PWKB=PWKB/sum(PWKB);

        end

        function [trajectory,times]=SSA(obj,maxTime,initialCellNumber)
        
        carryingCapacity=obj.carryingCapacity;    
        birthPropensity=obj.birthPropensity;
        birthPropensity(1)=[];
        deathPropensity=obj.deathPropensity;
        deathPropensity(end)=[];

        initialVectorSize=1e4;     
        times=nan(1,initialVectorSize+1);                                  %VECTOR OF TIMES OF OCURRENCE FOR EACH REACTION 
        trajectory=nan(1,initialVectorSize+1);                             %VECTOR OF STATES (OCCUPATION NUMBERS)         
        
        times(1)=0;                                                        %INITIAL TIME     
        trajectory(1)=initialCellNumber;                                   %INITIAL STATE  
        
        exitCondition=false;
        index=1;

        while exitCondition==false
    

            if index==numel(times)
               times=[times,nan(1,initialVectorSize+1)];
            end

            try
            rawPropensities=[birthPropensity(trajectory(index)),deathPropensity(trajectory(index))];%UPDATES PROPENSITIES
            propensityNormalisation=sum(rawPropensities);                  %UPDATES PROPENSITY NORMALISATION FACTOR     
            catch
            keyboard
            end

            draw=rand(1,2);                                                %PICKS TWO RANDOM NUMBERS UNIFORMLY DISTRIBUTED IN [0,1] 
            times(index+1)=times(index)+log(1/(1-draw(1)))/propensityNormalisation;%SETS TIME OF OCURRENCE FOR REACTION (DIRECT METHOD, COMES FROM INVERTING MOTE CARLO 
                                                                           %ALGOTIRHM)
            cumProbabilities=cumsum(rawPropensities/propensityNormalisation);%CHOOSES RANDOMLY (WEIGHTED BY NORMALISED PROPENSITIES) REACTION TO TAKE PLACE AND 
    
            I=find(cumProbabilities-draw(2)>=0,1,'first'); 
            trajectory(index+1)=(trajectory(index)+1)*(I==1)+(trajectory(index)-1)*(I==2);    

            if times(index+1)>=maxTime
                
                exitCondition=true;

            elseif trajectory(index+1)==0

                exitCondition=true;
                trajectory(index+2)=0;
                times(index+2)=maxTime;
            end

            index=index+1;

        end


        J=find(isnan(times),1,'first');
        
        if isempty(J)==0
            times(J:end)=[];
            trajectory(J:end)=[];
        end

        if sum(isnan(times))>0
            keyboard
        end

        end

        %Performs a clone size experiment for a birth-death process,
        %evaluating the clone size about the times within timeVector. The
        %numerical experiment is equivalent to labelling one cell and
        %follow the clonal evolution. To do so I consider the clone to be a
        %separate species (with initial condition of n=1) with the same
        %birth-death propensities, and simulate the two-species system
        %using Gillespie algorithm. 

        function cloneSizes=cloneSizeExperiment(obj,timeVector,initialCellNumber)
            
        %obj.birthPropensity(1)=[];
        obj.deathPropensity(end)=[];                                     
        
        cloneSizes=nan(1,numel(timeVector));
        %-----------------------INITIAL CONDITIONS-------------------------
     
        %initial clone number   
        cloneNumber=1; 
        %initial number of other cells
        cellNumber=initialCellNumber-1;
        %initial number of empty spaces
        emptyNumber=obj.carryingCapacity-initialCellNumber;      
        %initial time
        currentTime=0;
        %------------------------------------------------------------------
        
        %--------------ITERATION OF THE GILLESPIE ALGORITHM----------------
        %index for the time vector
        index=1; 
      
        while index<=numel(timeVector)

            %updates propensities    
            try rawPropensities=...
            [cloneNumber/(cloneNumber+cellNumber)*obj.birthPropensity(cloneNumber+cellNumber+1),...
            cloneNumber/(cloneNumber+cellNumber)*obj.deathPropensity(cloneNumber+cellNumber),...
            cellNumber/(cloneNumber+cellNumber)*obj.birthPropensity(cloneNumber+cellNumber+1),...
            cellNumber/(cloneNumber+cellNumber)*obj.deathPropensity(cloneNumber+cellNumber)];
            catch 
            keyboard
            end
            %updates normalisation factor
            propensityNormalisation=sum(rawPropensities);                     

            %picks two random numbers uniformly distributed in [0,1] 
            draw=rand(1,2);    
                
            %sets time in which a reaction takes place 
            currentTime=currentTime+log(1/(1-draw(1)))/propensityNormalisation;
            
            if currentTime > timeVector(index)
               cloneSizes(index)=cloneNumber;
               index=index+1;
            end
            
            %chooses which reaction is going to take place
            cumProbabilities=cumsum(rawPropensities/propensityNormalisation);
            I=find(cumProbabilities-draw(2)>0,1,'first'); 
            
            %updates cell numbers
            switch I
                case 1 %clone birth
                    cloneNumber=cloneNumber+1;
                    emptyNumber=emptyNumber-1;
                case 2 %clone death
                    cloneNumber=cloneNumber-1;
                    emptyNumber=emptyNumber+1;
                case 3 %other cell birth
                    cellNumber=cellNumber+1;
                    emptyNumber=emptyNumber-1;
                case 4 %other cell death
                    cellNumber=cellNumber-1;
                    emptyNumber=emptyNumber+1;      
            end
          
            %the following is just to speed up the simulation when n=0 is an
            %absorbing boundary. For other cases, comment these lines. 
            if cloneNumber==0 
                cloneSizes(index:end)=0;
                index=numel(timeVector)+1;
            end

        end        
       
        end

%          function cloneSizes=allCloneSizesvBDExperiment(obj,timeVector,initialCellNumber)
%              
%          obj.birthPropensity(1)=[];
%          obj.deathPropensity(end)=[];
%                  
%          numberOfClones=initialCellNumber;
%             
%          cloneSizes=nan(1,numel(timeVector));
%          %-----------------------INITIAL CONDITIONS------------------------
%      
%          %initial clone number   
%          cloneNumbers=ones(1,numberOfClones); 
%          %initial number of empty spaces
%          emptyNumber=obj.carryingCapacity-initialCellNumber;      
%          %initial time
%          currentTime=0;
%          %-----------------------------------------------------------------
%         
%         %--------------ITERATION OF THE GILLESPIE ALGORITHM----------------
%         %index for the time vector
%         index=1; 
%       
%         while index<=numel(timeVector)
%             
%             totalCellNumber=sum(cloneNumbers);
%              %updates propensities
%              for k=1:numberOfClones
%                  rawPropensities(k)=cloneNumbers(k)/totalCellNumber*obj.birthPropensity(totalCellNumber);
%                  rawPropensities(numberOfClones+k)=cloneNumbers(k)/totalCellNumber*obj.deathPropensity(totalCellNumber);
%              end
%             
%             %updates normalisation factor
%             propensityNormalisation=sum(rawPropensities);                     
% 
%             %picks two random numbers uniformly distributed in [0,1] 
%             draw=rand(1,2);    
%                 
%             %sets time in which a reaction takes place 
%             currentTime=currentTime+log(1/(1-draw(1)))/propensityNormalisation;
%             
%             if currentTime > timeVector(index)
%                cloneSizes(:,index)=cloneNumbers;
%                index=index+1;
% 
%             %the following is just to speed up the simulation when n=0 is an
%             %absorbing boundary. For other cases, comment these lines. 
%             elseif cloneNumber==0 
%                 cloneSizes(index:end)=0;
%                 index=numel(timeVector)+1;
%             end
% 
%             %chooses which reaction is going to take place
%             cumProbabilities=cumsum(rawPropensities/propensityNormalisation);
%             I=find(cumProbabilities-draw(2)>0,1,'first'); 
%             
%             %updates cell numbers
%             switch I
%                 case 1 %clone birth
%                     cloneNumber=cloneNumber+1;
%                     emptyNumber=emptyNumber-1;
%                 case 2 %clone death
%                     cloneNumber=cloneNumber-1;
%                     emptyNumber=emptyNumber+1;
%                 case 3 %other cell birth
%                     cellNumber=cellNumber+1;
%                     emptyNumber=emptyNumber-1;
%                 case 4 %other cell death
%                     cellNumber=cellNumber-1;
%                     emptyNumber=emptyNumber+1;      
%             end
%           
%         end        
%        
%          end

       
    end
end
