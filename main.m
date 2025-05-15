clear;
clc;
close all
dbstop if error

%Parameter setting
numMC=100;%Number of Monte Carlo Simulations
detection_prob=0.95;%Target missed detection rate
false_alarm_rate=60;%Clutter rate
drawP=0;%draw
GOSPA_plot=1;%GOSPA
PMBM_o=1;%GGIW_PMBM
PMBM_n=1;%PSM_MCGGIW_PMBM
isNewScenario=0;%generate new scenario

%Choose a scenario: 
%Scenario 1: Multiple irregular shapes maneuverable and intersection
%Scenario 2: Elliptical target maneuvering scenario
scenario = 1;

if scenario==1
    if isNewScenario
        star_create(numMC,detection_prob,false_alarm_rate,0);
    end
    modelparas1;
else
    if isNewScenario
        ellipse(numMC,detection_prob,false_alarm_rate,0);
    end
    modelparas2;
end

%Parameters used in GOSPA metric
c = 20;
p = 1;
%%%%%%%%%%%%%%%%%%%%%%  GGIW-PMBM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of time steps
K = model.K;
GOSPA_o = zeros(K,4,numMC);
trajectoryEstimates_o = cell(numMC,1);
simulation_time_o = zeros(length(Scenario.Z{1}),1);
if PMBM_o
    for t = 1:numMC
        simulation_time_this=0;
        Z = Scenario.Z{t};
        % Initialisation
        PPP.w = log(model.birth.w);
        PPP.GGIW = model.birth.GGIW;
        MBM.w = [];     % Global hypotheses weights
        MBM.track = {}; % Locl hypotheses trees
        MBM.table = []; % Global hypotheses look-up table
        
        estimates = cell(K,1);
        
        for k = 1:K
            
            tic
            %Print info
            pause(0);
         
            ['GGIW-PMBM    ',num2str(t),'    ',num2str(k)]
            %Update step
            [PPP,MBM] = updatePMBM(PPP,MBM,Z{k},k,model);
            
            %Extract estimates (both estimate of the current time and the estimate of the full trajectory)
            [estimates{k},trajectoryEstimates_o{t}{k}] = estimator(MBM,model);
            
            %Evaluate filtering performance using GOSPA 
            GOSPA_o(k,:,t) = GOSPAmetric(estimates{k},groundTruth{k},c,p);
            
            
            i0=size(estimates{k}.g,2);
            %Draw
            if drawP&&t==1
                figure(10+t);
                if scenario == 2
                    axis([100,600,100,600]);
                elseif scenario == 1
                    axis([100,700,100,700]);
                end

                hold on;
                for j = 1:i0
                    figure(10+t);
                    [x, y] = Sigmacircle_o(estimates{k}.x(1,j),estimates{k}.x(2,j),estimates{k}.X(:,:,j),2,3);
                end
            end
            drawnow;
            %Prediction Step
            if k < K
                [PPP,MBM] = predictPMBM(PPP,MBM,model); 
            end
            simulation_time_o(k) = simulation_time_o(k)+toc;
            simulation_time_this=simulation_time_this+toc;
        end
        disp(['time: ',num2str(simulation_time_this),' sec'])
    end
    GOSPA02_o = sum(GOSPA_o,3)/numMC;
    simulation_time_o = simulation_time_o/numMC;
    %GOSPA
    if GOSPA_plot
        for i=1:2
            figure(5);
            subplot(3,1,i)
            plot(1:K,GOSPA02_o(:,i),'-b');
            hold on;
        end
        subplot(3,1,3);
        plot(1:K,simulation_time_o(:),'-b');
        hold on;
    end
end


%%%%%%%%%%%%%%%%%%%%%%  PSM-MCGGIW-PMBM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GOSPA = zeros(K,4,numMC);
trajectoryEstimates = cell(numMC,1);
simulation_time = zeros(length(Scenario.Z{1}),1);

if PMBM_n
    for t = 1:numMC
        simulation_time_this=0;
        Z = Scenario.Z{t}; 
        PPP.w = log(model.birth.w);
        PPP.GGIW = model.birth.GGIW;
        MBM.w = [];     % Global hypotheses weights
        MBM.track = {}; % Locl hypotheses trees
        MBM.LineSegmentSet = {};%Multi-component edge expansion model
        MBM.table = []; % Global hypotheses look-up table
        
        estimates = cell(K,1);   
        for k = 1:K
            
            tic
            ['PSM-MCGGIW-PMBM    ',num2str(t),'    ',num2str(k)]
            %Update step
            [PPP,MBM] = updatePMBM_new(PPP,MBM,Z{k},k,model);
            
            %Extract estimates (both estimate of the current time and the estimate of the full trajectory)
            [estimates{k},trajectoryEstimates{t}{k}] = estimator(MBM,model);
            
            %Evaluate filtering performance using GOSPA
            GOSPA(k,:,t) = GOSPAmetric(estimates{k},groundTruth{k},c,p);
            
            
            i0=size(estimates{k}.g,2);
            if drawP&&t==1
                figure(10+t);
                if(~isempty(Z{k}))
                    plot(Z{k}(1,:),Z{k}(2,:),".black",'DisplayName','mes'); %draw measurements
                end
                hold on;
                for j = 1:i0
                    figure(10+t);
%                     [x, y] = Sigmacircle(estimates{k}.x(1,j),estimates{k}.x(2,j),estimates{k}.X(:,:,j),2,3);    %Overall extended plot
                    [~,idx] = max(MBM.w);
                    for l=1:length(MBM.LineSegmentSet{j}(MBM.table(idx,j)).LineSegment)
%                         plotSigmacircle(MBM.LineSegmentSet{j}(MBM.table(idx,j)).LineSegment{l}.Bern) %Extended component plot
                        plot([MBM.LineSegmentSet{j}(MBM.table(idx,j)).LineSegment{l}.start(1) MBM.LineSegmentSet{j}(MBM.table(idx,j)).LineSegment{l}.end(1)]...
                            ,[MBM.LineSegmentSet{j}(MBM.table(idx,j)).LineSegment{l}.start(2) MBM.LineSegmentSet{j}(MBM.table(idx,j)).LineSegment{l}.end(2)]...
                            ,'-','linewidth',1,'color','r','DisplayName','predictShape');%Marginal function plot
                    end
                end
                
            end
            drawnow;
            %Prediction Step
            if k < K
                
                [PPP,MBM] = predictPMBM_new(PPP,MBM,model);
            end
            simulation_time(k) = simulation_time(k)+toc;
            simulation_time_this=simulation_time_this+toc;
        end
        
        if drawP&&t==1
            if scenario == 2
                figure(10+t);
                axis equal
                axis([100,600,100,600]);
                
            elseif scenario == 1
                figure(10+t);
                axis equal
                axis([100,700,100,700]);
            end
        end
        disp(['time: ',num2str(simulation_time_this),' sec'])
    end
    GOSPA02 = sum(GOSPA,3)/numMC;
    simulation_time = simulation_time/numMC;
    if GOSPA_plot
        for i=1:2
            figure(5);
            subplot(3,1,i)
            plot(1:K,GOSPA02(:,i),'-r');
            hold on;
        end
        subplot(3,1,3);
        plot(1:K,simulation_time(:),'-r');
        hold on;
        subplot(3,1,1)
        title("GOSPA");
        subplot(3,1,2)
        title("LocationError");
        subplot(3,1,3)
        title("TimeCost");
    end
end



