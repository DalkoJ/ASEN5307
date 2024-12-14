%1D Kriging
%%
% Time Kriging
data_Vels=load('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\data_inSim.mat');
data_Cg=load("cg_sims.mat");
data_VelsDS=downsample(data_Vels.data_inSim,5); % Match the sim time.
data_CGs=data_Cg.cg_sim;
data_CGs=data_CGs(1:length(data_VelsDS),:);

casesWanted=[1 3 4 5 6 8 9 10 11 12 13 15 16]; %
casesTested=[2 7 14];

data_VelsDS_W=data_VelsDS(:,[1 casesWanted+1]);
data_CGs_W=data_CGs(:,[1 casesWanted+1]);

figure;
subplot(1,2,1)
for i=2:length(casesWanted)
    plot(data_VelsDS_W(1:201,1),data_VelsDS_W(1:201,i+1),'.'); hold on;
end
grid minor;
xlabel('Time [sec]');ylabel('Vel [m/sec]');
title('Direction Space Points')
subplot(1,2,2)
for i=2:length(casesWanted)
    plot(data_CGs_W(1:201,1),data_CGs_W(1:201,i+1),'.'); hold on;
end
grid minor;
ylabel('CG shift [m]');xlabel('Vel [m/sec]');
title('Output Space Points')
exportgraphics(figure(4),"avail_space_interpol.png")

statesRefs=[10 70 170 200];
[Y, Z] = meshgrid(-30:2:30, -1:.1:1); % Create a grid of y and z values

figure;
for i=statesRefs
    x=repelem(data_CGs_W(i,1),13)'; %Repeat element 13 times.
    % Define X as a constant value across the plane
    X = data_CGs_W(i,1) * ones(size(Y));  
    refPlane=surf(X, Y, Z,'LineStyle','none'); hold on; % Plot ref plane
    refPlane.FaceAlpha = 0.2;
    shading interp;
    scatter3(x,data_VelsDS_W(i,2:end)',...
        data_CGs_W(i,2:end)',20, data_CGs_W(i,2:end)', 'filled'...
        ,'MarkerEdgeColor', 'k');
end
zlim([-1 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
colorbar;
colormap('jet');
%title(["Original Design Space ","States: 10, 70, 170, 200","Times: 0.045, 0.345 0.845 0.995"])
title(["Original Design Space ","State: 170, 0.845 sec"])
exportgraphics(figure(2),"ex_st170_YZ_spline.png")
exportgraphics(figure(8),"ex_st170_YZ_differentCorrs.png")

%% compute variogram
Sbig=[];Ybig=[];
for i=170%:201%findSectionGroups
    Sloop=data_VelsDS_W(i,2:end)'; %From col 2 bc col 1 is time.
    Yloop=data_CGs_W(i,2:end)';

    Sbig=[Sbig; Sloop];
    Ybig=[Ybig; Yloop];
end

oneDvariogram(Sbig,Ybig) %Plot the Cloud variogram.
title("Cloud Variogram for State 170")% Variogram Estimator")

bookVariogram(Sbig,Ybig,'plotit',true)
grid minor;
xlabel('Averaged distance between observations')
ylabel('Averaged semivariance')
title("Experimental Variogram Estimator")

exportgraphics(figure(1),"ex_st170_CloudVariogram.png")

%%
interPolTime=1;
findSectionGroups=find(data_CGs_W(:,1)==interPolTime);
g1=2:1:51;
g2=52:1:101;
g3=102:1:151;
g4=152:1:201;
clear time_predict;
figure;
tlo=tiledlayout(10,5);
for i=170%:201%findSectionGroups
    Sloop=data_VelsDS_W(i,2:end)'; %From col 2 bc col 1 is time.
    Yloop=data_CGs_W(i,2:end)';
    %Sloop=xpp;
    %Yloop=ypp;
    [dmodel, perf]=dacefit(Sloop, Yloop, @regpoly1, @corrlin, .1, .01, .01);

    %X=[min(Sloop(:,1))*1.10:.01:max(Sloop(:,1))*1.10]';
    X=[-25:.01:25]';
    [YX ,MSE] = predictor(X, dmodel);
    time_predict{i}=[X, YX, MSE];
    %Find statistical parameters.
    %sigmaE=std(YX);
    upperT=YX+sqrt(MSE);%3*sigmaE;
    lowerT=YX-sqrt(MSE);%-3*sigmaE;
    
    nexttile;
    plot(X,upperT); hold on;
    plot(X,lowerT); hold on;
    %hold on; %Use this for the example
    %plot3(repelem(data_CGs_W(i,1),length(X))',X,YX,'b','Linewidth',1);
    %plot(X,YX+sigmaE,'--k'); hold on;
    %plot(X,YX-sigmaE,'--k'); hold on;
    plot(Sloop,Yloop,'* k','MarkerSize',5,'LineWidth',1.5); hold on;
    
    plot(X,YX,'-r','Linewidth',1.5);
    xlim([min(X) max(X)]);
    %xlim(-20 20]);
    %ylim([min(lowerT) max(upperT)]);
    ylim([-1 1]);
    %xlabel('Vel [m/sec]');ylabel('CG shift [m]');
    fill([X; flipud(X)], [upperT; flipud(lowerT)], 'b', ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Confidence interval');
    %title(["Prediction at t="+data_CGs(i,1)+" sec.", 'Corr. Model: Spline'])
    title("State: "+i+", "+data_CGs(i,1)+" sec.")
    %title(["State: "+i+", "+data_CGs(i,1)+" sec.",'Regr:0^{th}, Corr: Exp'])

    grid minor;
end
sgtitle("Regr: 1^{st}, Corr: Cubic Spline")
exportgraphics(figure(6),"g4_interpolFit.png")

save('predictedValsSpline.mat',"time_predict")
%plot3(repelem(data_CGs(i,1),16)',Sloop,Yloop,'x')

figure; 
surf(YX, X1,X2 ,'FaceColor','interp','edgecolor', 'k','LineStyle',':'); hold on;
colormap('turbo');
hold on;
plot3(repelem(data_CGs(i,1),16)',data_VelsDS(i,2:end)', data_CGs(i,2:end)','*k', 'MarkerSize',2); hold off;
%%

%%

gprMdl = fitrgp(Sloop,Yloop,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
ypred = resubPredict(gprMdl);

figure;
plot(Sloop,Yloop,'b.');
hold on;
plot(Sloop,Yloop,'--r','LineWidth',1.5);
xlabel('x');
ylabel('y');
legend('Data','GPR predictions');
hold off
    