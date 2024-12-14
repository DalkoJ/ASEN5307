%Develops the Kriging Interp.
%%
%data=load('data1.mat');
%S=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\TestS.mat'));
S=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\dace\S_space1.mat'));
S=S.S_space;
%S=S.tvar;
Y=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\dace\Y_space1.mat'));
%Y=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\dace\TestY.mat'));
Y=Y.Y_space;
%Y=Y.yvar;

figure;
scatter3(S(:,1),S(:,2), Y,20, Y, 'filled');
zlim([-1 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
colorbar;
colormap('jet');
title('Original Design Space')
exportgraphics(figure(1),"all_cg_points_3D.png")

%Compute Variogram
[X01,X02] = meshgrid(S(:,1));
[Y01,Y02] = meshgrid(S(:,2));
[Z01,Z02] = meshgrid(Y);

%Find linear distance
D = sqrt((X01 - X02).^2 + (Y01 - Y02).^2); 

% Experimental variogram.
G = 0.5*(Z01 - Z02).^2; 

indS(:,1) = 1:length(Y);
[C,R] = meshgrid(indS(:,1));
I = R > C;

%Plot variogram.
figure;
plot(D(I),G(I),'.b' )
xlabel('Lag distance')
ylabel('Semivariance')
title('Variogram Cloud')
grid minor;
exportgraphics(figure(2),"all_cg_points_cloud_Var.png")


%% Compute lag
D2 = D.*(diag(S(:,1)*NaN)+1);
lag = mean(min(D2));

hmd = max(D(:))/2;
max_lags = floor(hmd/lag);
LAGS = ceil(D/lag);

%Calculate classic variogram.
for i = 1 : max_lags
    SEL = (LAGS == i);
    DE(i) = mean(mean(D(SEL)));
    PN(i) = sum(sum(SEL == 1))/2;
    GE(i) = mean(mean(G(SEL)));
end

var_z = var(Y);
b = [0 max(DE)];
c = [var_z var_z];

figure;
plot(DE,GE,'.b','MarkerSize',15); hold on;
plot([0 max(DE)],[var_z var_z], '--k'); hold off;
xlim([0 max(DE)]);
ylim([0 1.1*max(GE)]) %10% tol.
xlabel('Averaged distance between seed points')
ylabel('Averaged semivariance')
title("Experimental Variogram Estimator")
grid minor;
exportgraphics(figure(3),"all_cg_points_exp_Var.png")

%% Separate data 
[posP,~]=find(Y>=0);
[posN,~]=find(Y<0);
S1=S(posP,:); S2=S(posN,:);
Y1=Y(posP); Y2=Y(posN);

figure;
scatter3(S1(:,1),S1(:,2), Y1,20, Y1, 'filled');
zlim([-1 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
colorbar;
colormap('jet');
title('Positive CG Design Space')
exportgraphics(figure(5),"pos_cg_points_3D.png")

figure;
scatter3(S2(:,1),S2(:,2), Y2,20, Y2, 'filled');
zlim([-1 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
colorbar;
colormap('jet');
title('Negative CG Design Space')
exportgraphics(figure(6),"neg_cg_points_3D.png")

%% Both Data
theta=[10 10];
lob=[.01 .01];
upb=[.1 .1];
[dmodel, perf]=dacefit(S, Y, @regpoly1, @corrgauss, theta, lob, upb);
%perf.nv %number of evaluations of the objective fnc.
% dmodel.theta %Theta-star
% dmodel.beta %Least Sq Estimate Beta-star
% dmodel.sigma2 %Process variance 

%% Find the dacefit for each Pos and Neg Z
[dmodelP, perfP]=dacefit(S1, Y1, @regpoly1, @corrgauss, theta, lob, upb);
[dmodelN, perfN]=dacefit(S2, Y2, @regpoly1, @corrgauss, theta, lob, upb);

%% Evaluate the predictor
%Create mapping grid.
sqMeshPts=250;
tolMap=.25; %[%]
%mapAreaLo=[min(S(:,1)) min(S(:,2))]*(1-tolMap);
%mapAreaHi=[max(S(:,1)) max(S(:,2))]*(1+tolMap);

mapAreaLo=[0 -30];
mapAreaHi=[5 30];

X = gridsamp([mapAreaLo; mapAreaHi], sqMeshPts);
[YX ,MSE] = predictor(X, dmodel);

%reshape the data.
X1 = reshape(X(:,1),sqMeshPts,sqMeshPts); 
X2 = reshape(X(:,2),sqMeshPts,sqMeshPts);
YX = reshape(YX, size(X1));

%% Evaluate pos predictor
%Create mapping grid.
[YXp ,MSEp] = predictor(X, dmodelP);
YXp = reshape(YXp, size(X1));%reshape the data.

[YXn ,MSEn] = predictor(X, dmodelN);
YXn = reshape(YXn, size(X1));%reshape the data.
%% Plot the Kriging.
figure; 
surf(X1, X2, YX,'FaceColor','interp','edgecolor', 'k','LineStyle',':'); hold on;
colormap('turbo'); 
plot3(S(:,1),S(:,2),Y,'*k', 'MarkerSize',2); hold off;
zlim([-1 1]); caxis([-1, 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
exportgraphics(figure(4),"krig_all.png")

figure; 
surf(X1, X2, YXp,'FaceColor','interp','edgecolor', 'k','LineStyle',':'); hold on;
surfP=gcf;
colormap('turbo'); 
plot3(S1(:,1),S1(:,2),Y1,'*k', 'MarkerSize',2); hold off;
zlim([-1 1]); caxis([-1, 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
title("Positive Component")

figure; 
surf(X1, X2, YXn,'FaceColor','interp','edgecolor', 'k','LineStyle',':'); hold on;
surfN=gcf;
colormap('turbo'); 
plot3(S2(:,1),S2(:,2),Y2,'*k', 'MarkerSize',2); hold off;
zlim([-1 1]); caxis([-1, 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('CG shift [m]');
grid minor;
title("Negative Component")

%% Overlay of the data
figure;
ax1 = axes;
surf(X1, X2, YXp,'FaceColor','interp','edgecolor', 'k','LineStyle',':'); hold on;
colormap(ax1, hot); % Apply "jet" colormap
alpha(0.7); % Adjust transparency
hold on;

% Create the second axes overlaying the first
ax2 = axes;
surf(X1, X2, YXn,'FaceColor','interp','edgecolor', 'k','LineStyle',':'); hold on;
colormap(ax2, cool); % Apply "hot" colormap
alpha(0.7); % Adjust transparency

% Link the axes and hide the top layer of ticks
% Make the second axes background transparent
ax2.Color = 'none';ax2.XColor = 'none';ax2.YColor = 'none';ax2.ZColor = 'none';

% Add labels
xlabel(ax1,'Time [sec]');
ylabel(ax1,'Vel [m/sec]');
zlabel(ax1,'CG shift [m]');
grid minor;

title(ax1, 'Positive and Negative Compontents');
colorbar(ax1, 'Position', [0.05, 0.1, 0.02, 0.8]); % Position for ax1 colorbar
colorbar(ax2, 'Position', [0.9, 0.1, 0.02, 0.8]); % Position for ax2 colorbar
%% Plot MSE
figure;
surf(X1, X2, reshape(MSE, size(X1)),'FaceColor','interp','edgecolor', 'w','LineStyle','none'); hold on;
colormap('turbo'); 
plot3(S(:,1),S(:,2),Y,'*w', 'MarkerSize',2); hold off;
zlim([-.001 1]); caxis([0, 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('Mean Sq Error [m]');
grid minor;
title("MSE Original Design Space")
exportgraphics(figure(8),"MSE_all.png")
%zlim([0 .01])
colorbar

figure;
surf(X1, X2, reshape(MSEp, size(X1)),'FaceColor','interp','edgecolor', 'w','LineStyle','none'); hold on;
colormap('turbo'); 
plot3(S(:,1),S(:,2),Y,'*k', 'MarkerSize',2); hold off;
zlim([-.001 1]); caxis([-.001, 1]);
xlabel('Time [sec]');ylabel('Vel [m/sec]');zlabel('Mean Sq Error');
grid minor;
%% Unwrapp the kriging maps
axesObjs = get(surfP, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%%get data 
XX1 = dataObjs.XData ;
YY1 = dataObjs.YData ;
ZZ1 = dataObjs.ZData ;

%Flatten the data.
XflatP = XX1(:);
YflatP = YY1(:);
ZflatP = ZZ1(:);
%% Finding testing points
St=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\S_test.mat'));
Stest=St.S_space;
Yt=load(fullfile('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\Y_test.mat'));
Ytest=Yt.Y_space;

for i=1:width(S)
    tpoint=[1.72691 8.79518];

    % Calculate Euclidean distances to the given point
    distances=sqrt((XflatP - tpoint(1,1)).^2 + (YflatP - tpoint(1,2)).^2);
    [~, closestIndex] = min(distances);% Find the index of the closest point

    % Extract the coordinates of the closest point
    closestPoint(i,:) = [Xflat(closestIndex), Yflat(closestIndex), Zflat(closestIndex)];
end

% Given point (x1, y1, z1)

% Display the closest point
disp('Closest point on the surface:');
disp(closestPoint);

% Highlight the closest point on the plot
plot3(closestPoint(1), closestPoint(2), closestPoint(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('Surface with Closest Point Highlighted');
xlabel('X');ylabel('Y');zlabel('Z');
hold off;