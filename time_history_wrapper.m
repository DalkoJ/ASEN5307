% Process all bearing data.
%% Choosing simulation cases
files_wanted_in=rmmissing(readlines(fullfile(pwd,filesep,'bearing_data',filesep,'cases_wanted.csv')));

for i=1:length(files_wanted_in)
    file_dir_loop=append(fullfile(pwd,filesep,'bearing_data'),filesep,files_wanted_in(i));
    accel_0i(:,i)=table2array(readtable(file_dir_loop));
end

%% Time series analysis
f_s=51.2*1000; %[Hz] Sampling frequency.
f_nyquist=f_s/2; %[Hz] Nyquis frequency.
t_0i=[0:1/f_s:15]'; %[sec] Time in
f_ds=1/.001; %[Hz] Downsampling frequency
t_in_max=5; %[sec] Time wanted for simulation.
range_wanted=1:find(t_0i==t_in_max); %Range of data wanted.
t_0=t_0i(range_wanted);
accel_0=accel_0i(range_wanted,:);

%min_peak_distance_idx = round(1000 / (f_s/length(accel_0Loop)));

figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Time hist    
    nexttile;
    plot(t_0,accel_0(:,i),'k'); hold on;
    xlim([0 t_in_max])
    grid minor;
    ylabel("Accel [G]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("Time Histories")
exportgraphics(figure(1),"raw_time_hist.png")

figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %downsampling
    [pxx_o,ff_o]=periodogram(accel_0(:,i),[],[],f_s); %Find frequency content
    
    nexttile;
    plot(ff_o,pxx_o,'k'); hold on;
    xlim([0 f_nyquist])
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
    
    %Find spectral peaks
    [pk{i}, lc{i}]=findpeaks(pxx_o,ff_o,'MinPeakDistance',500);  
    plot(lc{i},pk{i},'+r','MarkerSize',10) %Plot the spectral peaks.
end
title("Periodogram")
exportgraphics(figure(2),"raw_time_hist_period.png")


figure;
for i=1:length(files_wanted_in) 
    plot(lc{i},pk{i}); hold on;
end
grid minor;
ylabel("Power/Freq")
xlabel("Frequency [Hz]")
xlim([0 f_nyquist])
title("Spectral Peaks")
exportgraphics(figure(4),"raw_time_hist_spectral_peaks.png")

%% Detrend
plotThisSection='N';

for i=1:length(files_wanted_in)  
    
    H=[ones(length(accel_0(:,i)),1) t_0]; 
    %W=eye(length(accel_0(:,i))); %Weighting factors (all equal to 1).
    [U_b,S_b,V_b]=svd(H,0);
    s = diag(1./diag(S_b)); %Pseudoinverse of S.
    x_hatSV=V_b*s*U_b'*accel_0(:,i);
    y_fit=x_hatSV(1)+x_hatSV(2)*t_0;
    
%     a2=figure;
%     plot(t_0,accel_0(:,i),"x"); hold on;
%     plot(t_0,y_fit,'r','LineWidth',2);
%     grid minor;
%     legend("Observed","LSq Fit")
%     ylabel("Accel [G]")
%     xlabel("Time [sec]")
%     %exportgraphics(a2,"2a.png")
    
    %b: SE
    errorCov=inv(H'*H);
    x_hatInd=errorCov*H'*accel_0(:,i);
    sigma=diag(errorCov).^.5; %std error.
    accel_dt(:,i)=accel_0(:,i)-y_fit;
end 

figure;
yyaxis left
plot(skewness(accel_0),'xk','LineWidth',1); hold on;
plot(skewness(accel_dt),'o b','LineWidth',1); hold on;
ylabel("Skweness")
yline(0,'--')
yyaxis right
plot(kurtosis(accel_0),'+k','MarkerSize',10); hold on;
plot(kurtosis(accel_dt),'*r','LineWidth',1);
ylabel("Kurtosis")
xlabel("Cases")
yline(3,'--')
grid on;
title('Statistical Properties')
legend("Raw Sk.","Detr. Sk.","","Raw Ku.","Detr. Ku.",""...
    ,'Location','southoutside','Orientation','horizontal')
xticks([1:1:16])
exportgraphics(figure(7),"time_hist_Ku_Ske_dtr.png")

%% Low pass Filter design
figure;
for i=1:length(files_wanted_in)  
    plot(lc{i},pk{i}); hold on;
end
grid minor;
ylabel("Power/Freq")
xlabel("Frequency")
title("Spectral Peaks")

%n_th=4; %Trial
n_th=10; %Order of filter.
f_c= 5000; %Cutoff freq in [Hz]
[b,a]=butter(n_th,f_c/(f_s/2),'low');
[b1,a1]=cheby1(n_th,0.5,f_c/(f_s/2));
[h,w]=freqz(b,a);
[h1,w1]=freqz(b1,a1);

figure;
%hold on;
plot(w*0.5*f_s/pi, abs(h)); hold on;
plot(w1*0.5*f_s/pi, abs(h1)); hold on;
xline(f_c,':')
grid minor;
title("Filter Frequency Response")
xlim([0 f_nyquist])
legend('4th Order Butterworth','4th Order Chebychev1','10th Order Butterworth','10th Order Chebychev1')
xlabel('Frequency [Hz]')
exportgraphics(figure(3),"filter_comp.png")

%Filter the data.
clear accel_ft;
for i=1:length(files_wanted_in)  
    accel_ftLoop=filtfilt(b1,a1,accel_dt(:,i));
    accel_ft(:,i) = accel_ftLoop(nOrder/2 + 1:end - nOrder/2);
end
t_LF=t_0(nOrder/2 + 1:end - nOrder/2);

figure; %Plot their periodogram
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Find frequency content
    [pxx_o,ff_o]=periodogram(accel_ft(:,i),[],[],f_s); 
    [~, mI]=max(pxx_o);
    lowF_dominantF(i)=ff_o(mI);
    [~, lcLFLoop]=findpeaks(pxx_o,ff_o,'MinPeakDistance',50); 
    rotF=double(extractBetween(files_wanted_in(i),'N_','Hz'));
    if rotF==20
        indF=103; %Closest frequency to 40 Hz;
    elseif rotF==40
        indF=206; %Closest frequency to 40 Hz;
    else
        warning('Check rotF!')
    end
    [~, lcLFLoop2]=findpeaks(pxx_o,ff_o,'MinPeakHeight',pxx_o(indF)-.005); 
    lcLF{i}=lcLFLoop;
    lcLF2{i}=lcLFLoop2;

    nexttile;
    plot(ff_o,pxx_o,'k'); hold on;
    text(1200, max(pxx_o)*(1-.2), "f_{dom.}="+round(lowF_dominantF(i))+" Hz")%, 'FontSize', 12, 'Color', 'red');

    xlim([0 f_c/2])
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency [Hz]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
exportgraphics(figure(11),"filt_time_hist_period.png")

for i=1:length(files_wanted_in)
 LF_Table(i,:)=table(extractBefore(files_wanted_in(i),'.txt'),round((lcLF2{i}(1:5)')));
end
figure;
yyaxis left
plot(skewness(accel_0),'xk','LineWidth',1); hold on;
plot(skewness(accel_ft),'o b','LineWidth',1); hold on;
ylabel("Skweness")
yline(0,'--')
yyaxis right
plot(kurtosis(accel_0),'+k','MarkerSize',10); hold on;
plot(kurtosis(accel_ft),'*r','LineWidth',1);
ylabel("Kurtosis")
xlabel("Cases")
grid on;
title('Statistical Properties')
legend("Raw Sk.","Filt. Sk.","Raw Ku.","Filt. Ku."...
    ,'Location','southoutside','Orientation','horizontal')
xticks([1:1:16])
yline(3,'--')
exportgraphics(figure(6),"time_hist_Ku_Ske_filt.png")

%Peak Values
mean_accel_ft= mean(accel_ft);
std_noise_ft = std(accel_ft);
sigma_peaks_ft = mean_accel_ft' + 3 * std_noise_ft'; % 3-sigma threshold
skewness_ft=skewness(accel_ft);

for i=1:length(files_wanted_in)  
    [p0{1,i},~]=find(accel_ft(:,i)>sigma_peaks_ft(i));
    [p0{2,i},~]=find(accel_ft(:,i)<-sigma_peaks_ft(i));
end
%% Bandpass filter
f1=1.2*10^4;
f2=2.2*10^4;

nOrder= 70;
[b5,a5]=fir1(nOrder,[f1 f2]./(f_s/2),'bandpass');
[h5,w5]=freqz(b5,1,[],f_s);

figure; %Plot Filter
hold on;
xline(f1,"--"); hold on;
xline(f2,"--"); hold on;
plot(w5, abs(h5),'b','LineWidth',2); hold on;
grid minor;
ylabel("Magnitude")
xlabel("Frequency [Hz]")
title("Bandpass Filter Freq. Response")
legend(""+f1+"",""+f2+"","20th","30th","50th","70th"...
    ,'Location','southoutside','Orientation','horizontal')
exportgraphics(figure(12),"filter_comp_BP_different_ripples.png")

%Apply BP filter
for i=1:length(files_wanted_in) 
    accel_ft_HFLoop=filtfilt(b5,1,accel_dt(:,i));
    accel_ft_HF(:,i) = accel_ft_HFLoop(nOrder/2 + 1:end - nOrder/2);
end
t_HF=t_0(nOrder/2 + 1:end - nOrder/2);

figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    nexttile;
    plot(t_HF,accel_ft_HF(:,i),'k'); hold on;
    xlim([0 t_in_max])
    grid minor;
    ylabel("Accel [G]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("HF-Bandpass Time Histories")
exportgraphics(figure(14),"th_filter_comp_BP_end.png")

figure; %Plot their periodogram
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Find frequency content
    [pxx_o,ff_o]=periodogram(accel_ft_HF(:,i),[],[],f_s); 
    
    nexttile;
    plot(ff_o,pxx_o,'k'); hold on;
    xlim([0 f_s/2])
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
exportgraphics(figure(10),"th_filter_comp_BP_period.png")

%find peaks
for i=1:length(files_wanted_in)
    %Pos peaks every 1 second.
    [pkHF{i}, lcHF{i}]=findpeaks(accel_ft_HF(:,i),'MinPeakDistance',find(t_HF==1)); 
end

figure;
subplot(1,2,1)
for i=1:length(files_wanted_in)  
    plot(t_HF(lcHF{i}),pkHF{i},'+','Linewidth',2); hold on;
end
grid minor;
ylabel("Accel [G]")
xlabel("Time [sec]")
xlim([.95 5.05])
title("G Peaks")
subplot(1,2,2)
for i=1:length(files_wanted_in)  
    plot(t_HF(lcHF{i}),'+','Linewidth',2); hold on;
end
grid on;
xlabel("Peak Number")
ylabel("Time [sec]")
xlim([.95 5.05])
xticks([1:1:5])
title("Instance of Peak Gs")
exportgraphics(figure(8),"HF_gPeaks.png")

%% Extract HF pulses
accel_ft_pulses=[];
for i=1:length(files_wanted_in)
    %inM=find((t_HF-t_HF(1))==t_HF(lcHF{i}(1,1)));
    inM=lcHF{i}(1,1);
    accel_ft_pulses(:,i)=accel_ft_HF((inM-25):(inM+50),i);  
end
t_pulse=[0:1/f_s:(length(accel_ft_pulses)-1)/f_s]';
figure;
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    nexttile;
    plot(t_pulse,accel_ft_pulses(:,i),'b','linewidth',2); hold on;
    xlim([0 max(t_pulse)])
    grid minor;
    ylabel("Accel [G]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("Characteristic Pulse Shape")
exportgraphics(figure(17),"HF_pulses.png")

figure; %Plot their periodogram
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Find frequency content
    [pxx_HF,ff_HF]=periodogram(accel_ft_pulses(:,i),[],[],f_s);
    [~, mI]=max(pxx_HF);
    pulse_dominantF(i)=ff_HF(mI);

    nexttile;
    plot(ff_HF,pxx_HF,'b','linewidth',2); hold on;
    xlim([0 f_s/2])
    text(.3*10^4, max(pxx_HF)*(1-.2), "f_{dom.}="+pulse_dominantF(i)+" Hz")%, 'FontSize', 12, 'Color', 'red');

    %legend("f_{dominant} = "+pulse_dominantF(i)+" Hz", 'location','northwest')
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
exportgraphics(figure(16),"pulse_period.png")

figure; %Plot their periodogram
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Find frequency content
    [pxx_HF,ff_HF]=periodogram(accel_ft_pulses(:,i),[],[],f_s);
    [pxx_o,ff_o]=periodogram(accel_ft_HF(:,i),[],[],f_s); 
 
    nexttile;
    plot(ff_o,10*log10(pxx_o/f_s),'k'); hold on;
    plot(ff_HF,10*log10(pxx_HF/f_s),'b','linewidth',2);
    xlim([0 f_s/2])
    %text(.3*10^4, max(pxx_HF)*(1-.2), "f_{dom.}="+pulse_dominantF(i)+" Hz")%, 'FontSize', 12, 'Color', 'red');

    %legend("f_{dominant} = "+pulse_dominantF(i)+" Hz", 'location','northwest')
    grid minor;
    ylabel("dB")
    xlabel("Frequency [Hz]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
leg=legend("HF Data","Demodulated HF");
leg.Layout.Tile = 'south';
exportgraphics(figure(18),"demodulation.png")

%% Downsample LF of data for simulation
[p,q]=rat(f_ds/f_s);
t_ds=[min(t_LF):1/f_ds:max(t_LF)]';
plotThisSection='N';
%t_min_dist=.005; %This will alias some.
t_min_dist=.005;
accel_ds_forced=zeros(5000,16); flag1=0; accel_ds=[];
for i=1:length(files_wanted_in)      
    accel_ds(:,i)=resample(accel_ft(:,i),p,q);
    accel_ds_forced(1:end-1,i)=accel_ds(:,i);

    %MinPeakHeight with sigma_peaks does not work.
    [pk1{i}, lc1{i}]=findpeaks(accel_ft(:,i),f_s,'MinPeakDistance',t_min_dist); %Pos peaks.
    [pk2{i}, lc2{i}]=findpeaks(-accel_ft(:,i),f_s,'MinPeakDistance',t_min_dist); %Neg peaks.
    pk2{i}=-pk2{i}; %Invert the negative G value.
    
    if plotThisSection=='Y'
        aa=figure;
        plot(t_0,accel_ft(:,i),'k'); hold on;
        plot(t_ds,accel_ds(:,i),'r'); hold on;
        plot(lc1{i},pk1{i},'+y','MarkerSize',7,'LineWidth', 2); hold on;
        plot(lc2{i},pk2{i},'+y','MarkerSize',7,'LineWidth', 2);
        grid minor;
        ylabel("Accel [G]")
        xlabel("Time [sec]")
        title("Peak Spatial Distance D_p="+t_min_dist(i)+" [sec] - Peaks: +"+length(lc1)+"/-"+length(lc2)+"")
        %exportgraphics(aa,"tuning_dp_"+t_min_dist(i)+"_"+i+".png")  
    end

    %Find the nearest t_ds that corresponds to the peaks.
    [~, indMax] = min(abs(bsxfun(@minus, t_ds, lc1{i}.')), [], 1);
    [~, indMin] = min(abs(bsxfun(@minus, t_ds, lc2{i}.')), [], 1);
    
    if length(indMax)==length(indMin) %Should happen as an oscillatory manner.
        diff_index=indMax-indMin; %Ensures indeces don't overlap.
        flag1=0;
        if flag1==1
            figure;
            plot(diff_index);
            ylabel("Index Difference")
            xlabel("Peak Number")
            grid minor
            title("indMax-indMin Case:"+files_wanted_in(i)+"",'interpreter', 'none' )
        end 
        for j=1:length(diff_index)
            %Need to allocate values at different indeces.
            if diff_index(j)==0 %Nearest t_ds is the same.
                if lc1{1,i}(j,1)<lc2{1,i}(j,1)
                    %disp("Cond1: Max/min were found in the same index. Peak min/max "+j+" at "+lc2(j)+"/"+lc1(j)+"")
                    accel_ds_forced(indMax(j),i)=pk1{1,i}(j,1);
                    if indMax(j)==length(t_ds)
                        disp('HERE!')
                        continue
                    else
                        %Doesn't matter if referencing indMax or indMin.
                        %They should be the same at this instance.
                        accel_ds_forced(indMax(j)+1,i)=pk2{1,i}(j,1); %Assign the following value.
                    end
                    
                    if indMax(j)+1>length(t_ds) %If its at the end of tiem series, skip it.
                        continue %Don't do anything.
                    else
                        accel_ds_forced(indMax(j)+1,i)=pk1{1,i}(j,1); %Assign the following value.
                    end
                elseif lc1{1,i}(j,1)>lc2{1,i}(j,1)
                    %disp("Cond2: Max/min were found in the same index. Peak min/max "+j+" at "+lc2(j)+"/"+lc1(j)+"")
                    accel_ds_forced(indMin(j),i)=pk2{1,i}(j,1);
                    if indMax(j)+1>length(t_ds) %If its at the end of tiem series, skip it.
                        continue %Don't do anything.
                    else
                        accel_ds_forced(indMax(j)+1,i)=pk1{1,i}(j,1); %Assign the following value.                      
                    end
                else
                    %This should never happen but jic.
                    warning("lc1 and lc2 are the same!")
                end
            else 
                %Assign values.
                accel_ds_forced(indMax(j),i)=pk1{1,i}(j,1);
                accel_ds_forced(indMin(j),i)=pk2{1,i}(j,1);
            end
        end
    else
        warning("Lengths of superimposed indexes are not the same! Case:"+i+"")
        %Assign values.
        for j=1:length(indMax)
            accel_ds_forced(indMax(j),i)=pk1{1,i}(j,1);
        end
        for j=1:length(indMin)
            accel_ds_forced(indMin(j),i)=pk2{1,i}(j,1);
        end
    end
    %hold on;
    %plot(t_ds,accel_ds_forced,'b')
    %exportgraphics(aa,"input_accel"+t_min_dist+".png") 
end

%The forced interpolation may assign a value to the end of array.
accel_ds_forced=accel_ds_forced(1:end-1,:); %eliminate the last row.

figure; %Plot time histories downsampling.
%tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
     
    %nexttile;
    %plot(t_0,accel_dt(:,i),'k'); hold on;
    plot(t_LF,accel_ft(:,i),'b'); hold on;
    plot(t_ds,accel_ds(:,i),'g'); hold on;
    plot(t_ds,accel_ds_forced(:,i),'r'); hold on;
    xlim([0 t_in_max])
    grid minor;
    ylabel("Accel [G]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
sgtitle("Time Histories")
exportgraphics(figure(1),"downsampled_All.png")
legend("Filtered","Downsampled",'Imposed'...
    ,'Location','southoutside','Orientation','horizontal')

figure; %Plot their periodogram
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %Find frequency content
    [pxx_o,ff_o]=periodogram(accel_ds_forced(:,i),[],[],f_ds); 
    [pxx_o1,ff_o1]=periodogram(accel_ds(:,i),[],[],f_ds);
    [pxx_o0,ff_o0]=periodogram(accel_ft(:,i),[],[],f_s); 
    
    nexttile;
    loglog(ff_o0,10*log10(pxx_o0/f_s),'+b','linewidth',2); hold on
    loglog(ff_o1,10*log10(pxx_o1/f_ds),'g','linewidth',1); hold on
    loglog(ff_o,10*log10(pxx_o/f_ds),'r','linewidth',1);
    xlim([0 f_ds/2])
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end

figure; %Plot their periodogram
caseExample=16;

%Find frequency content
[pxx_o,ff_o]=periodogram(accel_ds_forced(:,caseExample),[],[],f_ds); 
[pxx_o1,ff_o1]=periodogram(accel_ds(:,caseExample),[],[],f_ds);
[pxx_o0,ff_o0]=periodogram(accel_ft(:,caseExample),[],[],f_s); 

figure;
subplot(3,1,1)
loglog(ff_o0,pxx_o0,'b','linewidth',2); hold on
xlim([1 f_ds/2])
grid minor;
ylabel("Power/Freq")
title("Filtered");
subplot(3,1,2)
loglog(ff_o1,pxx_o1,'g','linewidth',1); hold on
xlim([1 f_ds/2])
grid minor;
ylabel("Power/Freq")
title("Downsampled")
subplot(3,1,3)
loglog(ff_o,pxx_o/f_ds,'r','linewidth',1);
xlim([1 f_ds/2])
grid minor;
ylabel("Power/Freq")
xlabel("Frequency [Hz]")
title("Imposed")
sgtitle("Case: "+extractBefore(files_wanted_in(caseExample),'.txt')+"",'interpreter', 'none')
exportgraphics(figure(6),"downsampled_example_period.png")

%% Convert accel_ds_forced to velocity 
gravity_value=9.80665; %For earth.
vels=[];
for i=1:length(files_wanted_in)
    vels(:,i) = cumtrapz(accel_ds_forced(:,i)*gravity_value); % Velocity in Z
end 

figure; %Plot their periodogram
tlo=tiledlayout(4,4);
for i=16%:length(files_wanted_in)  
    %Find frequency content
    [pxx_o,ff_o]=periodogram(vels(:,i),[],[],f_ds); 
    
    %nexttile;
    loglog(ff_o,pxx_o,'k'); hold on;
    xlim([0 f_ds/2])
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
    xlim([0 100]);
end

exportgraphics(figure(3),"raw_vel_ex.png")


figure; %Plot time histories
%tlo=tiledlayout(4,4);
for i=16%:length(files_wanted_in)  
    %downsampling    
    nexttile;
    plot(t_ds,vels(:,i),'k'); hold on;
    xlim([0 t_in_max])
    grid minor;
    ylabel("Vel [m/sec]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
exportgraphics(figure(4),"raw_vel_exVel.png")


%% Fixing the velocity trends
f_cHF=100;
vels_ft=[];
for i=1:length(files_wanted_in) 
    if i==1  
        [vels_ft(:,i),filt_obj]=highpass(vels(:,i),f_cHF,f_ds,ImpulseResponse="iir",Steepness=0.95);
        [hV,wV]=freqz(filt_obj,1024,f_ds); %compute the filter once.
    else
        [vels_ft(:,i),~]=highpass(vels(:,i),f_cHF,f_ds,ImpulseResponse="iir",Steepness=0.95);
    end
end

figure;
plot(wV, abs(hV),'b','LineWidth',2); hold on;
grid minor;
ylabel("Magnitude")
xlabel("Frequency [Hz] ")
title("Filter Freq. Response")
xlim([0 500])
exportgraphics(figure(6),"vel_HP_filt.png")

figure; %Plot time histories
%tlo=tiledlayout(4,4);
for i=16%:length(files_wanted_in)  
    %downsampling    
    nexttile;
    plot(t_ds,vels_ft(:,i),'k'); hold on;
    xlim([0 t_in_max])
    grid minor;
    ylabel("Vel [m/sec]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
exportgraphics(figure(14),"vel_HP_filt_ex.png")

figure; %Plot their periodogram
%tlo=tiledlayout(4,4);
for i=16%:length(files_wanted_in)  
    %Find frequency content
    [pxx_o,ff_o]=periodogram(vels_ft(:,i),[],[],f_ds); 
    
    nexttile;
    plot(ff_o,pxx_o,'k'); hold on;
    xlim([0 f_ds/2])
    grid minor;
    ylabel("Power/Freq")
    xlabel("Frequency [Hz]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
    xlim([0 500]);
end
exportgraphics(figure(13),"vel_period_filt_ex.png")


vels_dt=[];
for i=1:length(files_wanted_in)      
    H=[ones(length(vels_ft(:,i)),1) t_ds]; 
    [U_b,S_b,V_b]=svd(H,0);
    s = diag(1./diag(S_b)); %Pseudoinverse of S.
    x_hatSV=V_b*s*U_b'*vels_ft(:,i);
    y_fit=x_hatSV(1)+x_hatSV(2)*t_ds;
    
    %b: SE
    errorCov=inv(H'*H);
    x_hatInd=errorCov*H'*vels_ft(:,i);
    sigma=diag(errorCov).^.5; %std error.
    vels_dt(:,i)=vels_ft(:,i)-y_fit;    
end 

%Reduced Peak Veloecities beyond 3sigma Values
mean_vel= mean(vels_dt);
std_vel = std(vels_dt);
sigma_vel = mean_vel' + 3 * std_vel'; % 3-sigma threshold
skewness_ft=skewness(vels_dt);
for i=1:length(files_wanted_in) 
    vels_dtLoop=vels_dt(:,i); %Assign new variable to preserve old vel.
    vels_dtLoop(find(vels_dtLoop>sigma_vel(i)*1.1))=sigma_vel(i); %+ve peaks
    vels_dtLoop(find(vels_dtLoop<-sigma_vel(i)))=-sigma_vel(i); %-ve peaks
    vels_dt_s(:,i)=vels_dtLoop;
end

figure; %Plot time histories
tlo=tiledlayout(4,4);
for i=1:length(files_wanted_in)  
    %downsampling    
    nexttile;
    %plot(t_ds,vels(:,i),'k'); hold on;
    %plot(t_ds,vels_ft(:,i),'b'); hold on;
    plot(t_ds,vels_dt(:,i),'+r'); hold on;
    plot(t_ds,vels_dt_s(:,i),'m'); hold on;
    xlim([0 t_in_max])
    yline(sigma_vel(i),'--'); hold on;
    yline(-sigma_vel(i),'--');
    grid minor;
    ylabel("Vel [m/sec]")
    xlabel("Time [sec]")
    title("Case: "+files_wanted_in(i)+"",'interpreter', 'none' )
end
leg=legend("Filtered",'No 3\sigma','Orientation','Horizontal');
leg.Layout.Tile = 'south';
exportgraphics(figure(16),"final_velocity.png")

%% Wrap accel_in to .k file for simulation
formatK=rmmissing(readlines("accels_format.k"));
filename = 'accels.k';
mkDirFlag='N';
for i=1:length(files_wanted_in) 
    clear k_file
    caseName=extractBefore(files_wanted_in(i),'.txt');
    if mkDirFlag=='Y'
        mkdir(""+caseName+"")
    end

    fileIDScrap = fopen('accels_scrap.k','wt');
    for j=1:length(t_ds)
        fprintf(fileIDScrap,'%-19d %-19d\n',t_ds(j), accel_ds_forced(j,i))
    %fprintf(fileIDScrap, '%-9s %-9s\n', t_ds, accel_ds_forced(:,i))
    end
    fclose(fileIDScrap)
    
    %Read printed values and add "END"
    accelScrap=readlines("accels_scrap.k");
    accelScrap(end)="*END";
    
     % Merge with baseline boilerplate.
    k_file=[formatK;accelScrap];
    
    fileID = fopen(fullfile(pwd,filesep,caseName,'accels.k'),'wt');
    fprintf(fileID,'%s\n',k_file)
    fclose(fileID)
end

%% Wrap HF vels_dt to .k file for simulation
formatK=rmmissing(readlines("vels_format.k"));
filename = 'vels.k';
mkDirFlag='Y';
for i=1:length(files_wanted_in) 
    clear k_file
    caseName=strcat(extractBefore(files_wanted_in(i),'.txt'),"_LF");
    if mkDirFlag=='Y'
        mkdir(""+caseName+"")
    end

    fileIDScrap = fopen('format_scrap.k','wt');
    for j=1:length(t_ds)
        fprintf(fileIDScrap,'%-19d %-19d\n',t_ds(j), vels_dt_s(j,i))
    end
    fclose(fileIDScrap)
    
    %Read printed values and add "END"
    accelScrap=readlines("format_scrap.k");
    accelScrap(end)="*END";
    
     % Merge with baseline boilerplate.
    k_file=[formatK;accelScrap];
    
    fileID = fopen(fullfile(pwd,filesep,caseName,'vels.k'),'wt');
    fprintf(fileID,'%s\n',k_file)
    fclose(fileID)
end

data_inSim=[t_ds, vels_dt_s];
