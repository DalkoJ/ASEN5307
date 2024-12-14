%Tests kriging interpo
%%
mainFolder='C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\';
files_wanted_in=rmmissing(readlines(fullfile(mainFolder,filesep,'bearing_data',filesep,'cases_wanted.csv')));

%Estimate CGS
pVals1=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\dace\predictedValsAll.mat'));
pVals2=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\dace\predictedVals.mat'));
pVals3=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\dace\predictedValsSpline.mat'));

%pVals=pVals.time_predict; %Each col is a time step interpol.
pVals(1,:)=pVals1.pValsOld; %All Gaussian
pVals(2,:)=pVals2.time_predict; %Gaussian
pVals(3,:)=pVals3.time_predict; %Spline

predictorNames=["All Gaussian", "Gaussian", "Cubic Spline"];
%Input pars 
data_Vels=load('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\data_inSim.mat');
data_VelsDS=downsample(data_Vels.data_inSim,5); % Match the sim time.

%Simulated CGs
data_Cg=load("cg_sims.mat");
data_CGs=data_Cg.cg_sim; %Colums are each test case vals.
data_CGs=data_CGs(1:length(data_VelsDS),:);

%Assemble input vars.
casesTested=[2 7 14];
data_VelsDS_T=data_VelsDS(:,[1 casesTested+1]);
data_CGs_T=data_CGs(:,[1 casesTested+1]); %Used for verification.

finishTimeI=find(data_VelsDS_T(:,1)==1); %==1 for 1 second.
%% 
clear cg_rtVals MSE_rtVals cg_rtValsAll MSE_rtValsAll
for k=1:3
    for j=1:length(casesTested)
        for i=1:finishTimeI            
            preValsLoop=pVals{k,i};
            [~, closestIndex] = min(abs(data_VelsDS_T(i,j+1)-preValsLoop(:,1)));
            
            if closestIndex==1 || closestIndex==length(preValsLoop)
                warning("Value is outside the Kriging interpolation.")
                cg_rtVals(i)=NaN;
                MSE_rtVals(i)=NaN;
            else
                %Retrieve col 2 and 3 (CG-est and MSE)
                cg_rtVals(i)=preValsLoop(closestIndex,2); 
                MSE_rtVals(i)=preValsLoop(closestIndex,3);
            end
        end
        cg_rtValsLoop(:,j)=cg_rtVals;
        MSE_rtValsLoop(:,j)=MSE_rtVals;
    end
    cg_rtValsAll{k}=cg_rtValsLoop;
    MSE_rtValsAll{k}=MSE_rtValsLoop;
end
%%
figure;
tiledlayout(3,3)
for k=1:3
    for i=1:length(casesTested)
        nexttile
        cg_rtVals=cg_rtValsAll{k};
        plot(data_CGs_T(1:finishTimeI,1),data_CGs_T(1:finishTimeI,i+1)...
            ,'k', 'LineWidth',2);hold on;
        plot(data_VelsDS_T(1:finishTimeI,1),cg_rtVals(:,i),'b', 'LineWidth',.75); 
        
        title(["Cases"+extractBefore(files_wanted_in(casesTested(i)),'.txt')+""...
            ,"Corr Model: "+predictorNames(k)+""],'Interpreter','none')
        grid minor
        xlabel('Time [sec]');ylabel('CG shift [m]');
        ylim([-1 1])
    end
end        
leg=legend("Simulated","Interpolated",'Orientation','horizontal');
leg.Layout.Tile = 'south';

exportgraphics(figure(10),"final_predictor_comparison.png")
%% Residuals
figure;
tiledlayout(3,3)
for k=1:3
    for i=1:length(casesTested)
        nexttile
        cg_rtVals=cg_rtValsAll{k};
        residualsLoop=cg_rtVals(:,i)-data_CGs_T(1:finishTimeI,i+1);
        residualsAll{k,i}=residualsLoop;
        normplot(residualsLoop)
        grid minor;
         
        title(["Cases"+extractBefore(files_wanted_in(casesTested(i)),'.txt')+""...
            ,"Corr Model: "+predictorNames(k)+"","Residuals Ku:"+kurtosis(residualsLoop)+""],'Interpreter','none')
        grid minor
        xlabel('Probability');xlabel('CG shift Residuals[m]');
        %ylim([-1 1])
    end
end 
exportgraphics(figure(2),"final_predictor_residuals_nomplot.png")


t_in=data_VelsDS_T(1:finishTimeI,1);
figure;
tiledlayout(3,3)
for k=1:3
    for i=1:length(casesTested)
        
        t_inLoop=t_in;
        residualsLoop=residualsAll{k,i};
        flagN=isnan(residualsLoop(:,1));
        residualsLoop(flagN) = [];
        t_inLoop(flagN) = [];
        nexttile
        plot(t_inLoop,residualsLoop,".b"); hold on;
        for i=1:3
            yline(i*std(residualsLoop),'--k'); hold on;
            yline(-i*std(residualsLoop),'--k'); hold on;
        end
        ylabel("Residuals (Sim.-Interp.) [m]")
        xlabel("Time [sec]")
        title(["Cases"+extractBefore(files_wanted_in(casesTested(i)),'.txt')+""...
            ,"Corr Model: "+predictorNames(k)+"",...
            "Residuals Ku:"+kurtosis(residualsLoop)+""],'Interpreter','none')
        grid minor
    end
end
exportgraphics(figure(9),"final_predictor_residuals_th.png")



