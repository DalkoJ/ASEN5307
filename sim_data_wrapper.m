%Process all simulation data.
mainFolder='C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\';
files_wanted_in=rmmissing(readlines(fullfile(mainFolder,filesep,'bearing_data',filesep,'cases_wanted.csv')));
bearing_geom=readtable(fullfile(mainFolder,filesep,'bearing_geom.txt'));
gapIDval=string(cell2mat(bearing_geom{:,1}));

simFolder=fullfile(mainFolder,filesep,'simData5sec');
for i=1:length(files_wanted_in)
    caseName=extractBefore(files_wanted_in(i),'.txt');
    nFile1=strcat(caseName,"_LF_float.txt");
    nFile2=strcat(caseName,"_LF_cgZ.txt");
    nFile3=strcat(caseName,"_LF_pd.txt");

    %Extract floating level sim data.
    float_zLoop=table2array(readtable(append(fullfile(simFolder,...
        filesep,nFile1))));
    cg_zLoop=table2array(readtable(append(fullfile(simFolder,...
        filesep,nFile2))));
    pd_Loop=table2array(readtable(append(fullfile(simFolder,...
        filesep,nFile3))));
    %find the mean Pressure and convert to kPa
    pd_mean(:,i)=mean(pd_Loop(:,2:end),2)*10^(-3); 

    %Extract cg level sim data.
    float_z(:,i)=float_zLoop(:,2);
    cg_x(:,i)=cg_zLoop(:,2);
    cg_y(:,i)=cg_zLoop(:,3);
    cg_z(:,i)=cg_zLoop(:,4);
 
    %Extract time once for each sim data.
    if i==1
        t_float=float_zLoop(:,1);
        t_cg=cg_zLoop(:,1);
        t_pd=pd_Loop(:,1);
    end
end

%Plotting data in
figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Time hist    
    nexttile;
    plot(t_cg,cg_x(:,i),'k'); hold on;
    plot(t_cg,cg_y(:,i),'g'); hold on;
    plot(t_cg,cg_z(:,i),'b','LineWidth',1.5); hold on;
    xlim([0 max(t_cg)])
    grid minor;
    ylabel("CG shift [m]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("CG Time Histories")
leg=legend("X-Axis",'Y-Axis','Z-Axis','Orientation','Horizontal');
leg.Layout.Tile = 'south';
exportgraphics(figure(10),"all_raw_cg.png")

figure; %Plot floater
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Time hist    
    nexttile;
    plot(t_float,float_z(:,i),'b'); hold on;
    yline(.25,'--k','lineWidth',1)
    xlim([0 max(t_float)])
    grid minor;
    ylabel("Floating Point [m]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("Floating Point Time Histories")

figure; %Plot Nodal Pressure
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Time hist    
    nexttile;
    plot(t_pd,pd_mean(:,i),'b'); hold on;
    yline(.25,'--k','lineWidth',1)
    xlim([0 max(t_float)])
    grid minor;
    ylabel("Mean Pressure [kPa]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("Mean Pressure Time Histories")
%% Post process data.
cg_rss=[];
for i=1:length(files_wanted_in)
    cg_rss(:,i)=(cg_x(:,i).^2+cg_y(:,i).^2+cg_z(:,i).^2).^0.5; %Find the RSS delta-CG
    vel_cg(:,i)=diff(cg_rss(:,i))./diff(t_cg);
    %acel_cg(:,i)=diff(vel_cg(:,i))./diff(t_cg);
end

% Highpass Filter the data to remove oscillatory trends.
f_s1=1/((t_cg(3)-t_cg(2)));
f_cHF=1;
cg_rss_ft=[];
for i=1:length(files_wanted_in) 
    if i==1  
        [cg_rss_ft(:,i),filt_obj]=highpass(cg_rss(:,i),f_cHF,f_s1,ImpulseResponse="iir",Steepness=0.95);
        [hV,wV]=freqz(filt_obj,1024,f_s1);
    else
        [cg_rss_ft(:,i),~]=highpass(cg_rss(:,i),f_cHF,f_s1,ImpulseResponse="iir",Steepness=0.95);
    end
end

figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %downsampling
    [pxx_o,ff_o]=periodogram(cg_rss_ft(:,i),[],[],f_s1); %Find frequency content
    [~,mI(i)]=max(pxx_o);

    nexttile;
    plot(ff_o,pxx_o,'k'); hold on;
    xlim([0 5]) %5 was chosen after seeing up to the Nyquist.
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency [Hz]")
    title(["Case: "+files_wanted_in(i)+"","Dom. F:"+ff_o(mI(i))+""],'interpreter', 'none' )
end
title("CG Periodogram")
exportgraphics(figure(9),"all_cg_period.png")


figure;
plot(wV, abs(hV),'b','LineWidth',2); hold on;
grid minor;
ylabel("Magnitude")
xlabel("Frequency [Hz] ")
title("Filter Freq. Response")
xlim([0 3])
exportgraphics(figure(9),"cg_HP_filt.png")


%Detrend CG
cg_dt=[];
for i=1:length(files_wanted_in)      
    H=[ones(length(cg_rss_ft(:,i)),1) t_cg]; 
    [U_b,S_b,V_b]=svd(H,0);
    s = diag(1./diag(S_b)); %Pseudoinverse of S.
    x_hatSV=V_b*s*U_b'*cg_rss_ft(:,i);
    y_fit=x_hatSV(1)+x_hatSV(2)*t_cg;
    
    errorCov=inv(H'*H);
    x_hatInd=errorCov*H'*cg_rss_ft(:,i);
    sigma=diag(errorCov).^.5; %std error.
    cg_dt(:,i)=cg_rss_ft(:,i)-y_fit; 
end 
cg_sim=[t_cg, cg_dt];
save("cg_sims.mat",'cg_sim')

figure; %Plot cg
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)   
    nexttile;
    %plot(t_cg,cg_rss(:,i),'b'); hold on;
    plot(t_cg,cg_dt(:,i),'k'); hold on;
    xlim([0 max(t_cg)])
    grid minor;
    ylabel("Position [m]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("CG RSS Time Histories")
exportgraphics(figure(9),"cg_all_filt.png")

figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %downsampling
    [pxx_o,ff_o]=periodogram(cg_dt(:,i),[],[],f_s1); %Find frequency content
    [~,mI(i)]=max(pxx_o);

    nexttile;
    plot(ff_o,pxx_o,'k'); hold on;
    xlim([0 5]) %5 was chosen after seeing up to the Nyquist.
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency [Hz]")
    title(["Case: "+files_wanted_in(i)+"","Dom. F:"+ff_o(mI(i))+""],'interpreter', 'none' )
end
sgtitle("CG Periodogram")
exportgraphics(figure(12),"cg_all_per_filt.png")

figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %downsampling
    figure(i);
    cwt(cg_dt(:,i),f_s1); %Find frequency content
    
    grid minor;
    xlabel("Time [sec]")
    ylabel("Frequency [Hz]")
    title(["Case: "+files_wanted_in(i)+"","Dom. F:"+ff_o(mI(i))+""],'interpreter', 'none' )
    exportgraphics(figure(i),"cwt_"+i+"_filt.png")
end
sgtitle("CG Periodogram")

%Integrate to velocity
for i=1:length(files_wanted_in)
    vel_cg(:,i)=diff(cg_dt(:,i))./diff(t_cg);
end

figure; %Plot Vel RSS
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Time hist    
    nexttile;
    plot(t_cg(1:end-1),vel_cg(:,i),'b'); hold on;
    %yline(.25,'--k','lineWidth',1)
    xlim([0 t_cg(end-1)])
    grid minor;
    ylabel("Floating Point [m]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("Velocity of CG_{RSS} Time Histories")

%% Obtain S design Space
load('C:\Users\dalko\OneDrive\Desktop\ASEN 5307\final_project\data_inSim.mat');

data_inSimDS=downsample(data_inSim,5); % Match the sim time.
t_min_dist=.05;
f_s2=(data_inSimDS(3,1)-data_inSimDS(2,1))^-1;
pkp=[]; lkp=[]; cg_p=[];
pkn=[]; lkn=[]; cg_n=[];
for i=1:length(files_wanted_in)  
    %MinPeakHeight with sigma_peaks does not work.
    [pkpLoop, lkpLoop]=findpeaks(data_inSimDS(:,i+1),'MinPeakDistance',10); %Pos peaks.
    [pknLoop, lknLoop]=findpeaks(-data_inSimDS(:,i+1),'MinPeakDistance',10); %Neg peaks.
    pknLoop=-pknLoop; %Invert the negative G value. 
    pkp=[pkp; pkpLoop];
    pkn=[pkn;pknLoop];
    lkp=[lkp; lkpLoop];
    lkn=[lkn; lknLoop];

    %Find cg respective values
    cg_pLoop=cg_dt(lkpLoop,i);
    cg_p=[cg_p; cg_pLoop];
    cg_nLoop=cg_dt(lknLoop,i); 
    cg_n=[cg_n; cg_nLoop];
end

pkf=[pkp; pkn];
lkf=data_inSimDS([lkp; lkn],1);
cgf=[cg_p; cg_n];

figure;
%plot3(lkf,pkf,cgf,'.'); hold on;
plot(lkf,pkf,'.'); hold on;
yline(std(pkf,0,'all'),'--k'); hold on;
yline(-std(pkf,0,'all'),'--k'); hold on;
yline(2*std(pkf,0,'all'),'--k'); hold on;
yline(-2*std(pkf,0,'all'),'--k'); hold off;

d_hin=cgf;
bRangeData=linspace(min(d_hin),max(d_hin));
bDist=makedist("Normal",mean(d_hin,'all'),std(d_hin,0,'all'));
bPdf=pdf(bDist,bRangeData); 
nBins1=round(length(d_hin)^0.5);
[b1,b2]=chi2gof(d_hin,'alpha',0.05);
if b1==1
    disp("   - Distribution DOES NOT follows a Normal distribution")
elseif b1==0
    disp("   - Distribution DOES follows a Normal distribution")
else
    warning("Check this chi2gof.m call")
end
figure;
bHist=histogram(d_hin,nBins1,'Normalization','pdf'); hold on;
plot(bRangeData,bPdf,'b','LineWidth',2); hold on;
xline(std(d_hin,0,'all'),'--k'); hold on;
xline(-std(d_hin,0,'all'),'--k'); hold on;
xlabel("Data"); ylabel("Probability [100^{-1}]")
title(["cg_f Distribution (Bins="+nBins1+")","Ku="+kurtosis(d_hin)+""])
grid minor;

%% Shaping per cgF 1 Used in finding peak values. NOT USED IN FINAL RESULTS
pkn2=[]; lkn2=[]; lkp2=[]; cg_n2=[]; cg_p2=[]; pkp2=[]; 
thresholVal=std(cgf,0,'all')/1.25;
thresholVal=std(cgf,0,'all')/1.25;
d_inT=cgf;
%MinPeakHeight with sigma_peaks does not work.
[pos,~]=find(d_inT>-thresholVal & d_inT<thresholVal);
[cg_p2, lkp2]=findpeaks(d_inT(pos),'MinPeakDistance',4); %Pos peaks.
[cg_n2, lkn2]=findpeaks(-d_inT(pos),'MinPeakDistance',4); %Neg peaks.
cg_n2=-cg_n2; %Invert the negative G value. 

pkf_Loopt1=pkf([lkp2;lkn2]);

%Assemble remaining values.
cgf_remain=cgf(setdiff(1:length(cgf), pos));
pkf_remain=pkf(setdiff(1:length(pkf), pos));
lkf_remain=lkf(setdiff(1:length(lkf), pos));

cgf_t1=[cgf_remain;cg_p2;cg_n2];
pkf_t1=[pkf_remain;pkf_Loopt1];
lkf_t1=[lkf_remain;lkf([lkp2;lkn2])];

d_hin=cgf_t1;
bRangeData=linspace(min(d_hin),max(d_hin));
bDist=makedist("Normal",mean(d_hin,'all'),std(d_hin,0,'all'));
bPdf=pdf(bDist,bRangeData); 
nBins1=round(length(d_hin)^0.5);
[b1,b2]=chi2gof(d_hin,'alpha',0.01);
if b1==1
    disp("   - Distribution DOES NOT follow a Normal distribution")
elseif b1==0
    disp("   - Distribution DOES follow a Normal distribution")
else
    warning("Check this chi2gof.m call")
end
figure;
bHist=histogram(d_hin,nBins1,'Normalization','pdf'); hold on;
plot(bRangeData,bPdf,'b','LineWidth',2); hold on;
xline(std(d_hin,0,'all'),'--k'); hold on;
xline(-std(d_hin,0,'all'),'--k'); hold on;
xlabel("Data"); ylabel("Probability [100^{-1}]")
title(["cg_f Distribution (Bins="+nBins1+")","Ku="+kurtosis(d_hin)+""])
grid minor;

figure;
%plot3(lkf,pkf,cgf,'.'); hold on;
plot(lkf_t1,pkf_t1,'.'); hold on;
yline(std(pkf_t1,0,'all'),'--k'); hold on;
yline(-std(pkf_t1,0,'all'),'--k'); hold on;
yline(2*std(pkf_t1,0,'all'),'--k'); hold on;
yline(-2*std(pkf_t1,0,'all'),'--k'); hold off;
%exportgraphics(figure(5),"cgtf_dist_t6_1s.png")


%% Preassemble S and Y
files_wanted_k=rmmissing(readlines(fullfile(mainFolder,filesep,'cases_wantedMod.csv')));
files_wanted_f=[];files_tested_f=[];
for i=1:length(files_wanted_k)
    %Skip files that are used for testing the model.
    if contains(files_wanted_k(i),'$')
        files_tested_floop=extractAfter(files_wanted_k(i),'$');  
        files_tested_f=[files_tested_f; files_tested_floop];
    else
        files_wanted_fLoop=files_wanted_k(i);
        files_wanted_f=[files_wanted_f; files_wanted_fLoop];
    end
end

%Find indices that correspond to the simData location.
respectiveWantedFiles=[]; respectiveTestedFiles=[];
for i=1:length(files_wanted_f)
    respectiveWantedFiles(i,1)=find(files_wanted_f(i)==files_wanted_in);
end
for i=1:length(files_tested_f)
    respectiveTestedFiles(i,1)=find(files_tested_f(i)==files_wanted_in);
end
%% Assemble S and Y
flagAssembly='W'; %INPUT: 'T': tested

if flagAssembly=='W'
    respectiveFiles=respectiveWantedFiles;
    files_neutral=files_wanted_f;
elseif flagAssembly=='T'
    respectiveFiles=respectiveTestedFiles;
    files_neutral=files_tested_f;
else
    warning("Flag must 'T' or 'W'.")
end

Y_space=cgf_t1;
Y_space=[]; 
for i=1:length(files_neutral) 
    %Y_space(i,1)=mean(float_z(:,respectiveFiles(i)));
    %Y_space(i,2)=max(vel_cg(:,respectiveFiles(i)));
    %Y_space(i,1)=cg_dt(cgf,respectiveFiles(i));
end
S_space=[lkf_t1, pkf_t1];
S_space=[];
for i=1:length(files_neutral)
    gapIdLoop=extractBefore(files_neutral(i),'_');

    S_space(i,1:2)=bearing_geom{find(gapIdLoop==gapIDval),2:3};
    S_space(i,3)=double(extractBetween(files_neutral(i),'_','N_'));
    S_space(i,4)=double(extractBetween(files_neutral(i),'N_','Hz'));
end

if flagAssembly=='W'
    save('S_space1.mat','S_space')
    save('Y_space1.mat','Y_space')
elseif flagAssembly=='T'
    save('S_test.mat','S_space')
    save('Y_test.mat','Y_space')
end
