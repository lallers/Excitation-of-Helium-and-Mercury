%% Helium Excitation
close all;clc;clear all
chr = strcat(mfilename(),'.m');
path = which(mfilename());path = strrep(path,chr,'');

he_r1_v248 = load(strcat(path,'data\DAY_2_RUN1_V248.txt'));
he_r2_v248 = load(strcat(path,'data\DAY_2_RUN2_V248.txt'));
he_r3_v248 = load(strcat(path,'data\DAY_2_RUN3_V248.txt'));
he_r4_v248 = load(strcat(path,'data\DAY_2_RUN4_V248.txt'));
test = load(strcat(path,'data\run3_4peaks.txt'));
test2 = load(strcat(path,'data\run6_revesedbias.txt'));
he_v284.vol = [he_r1_v248(:,1), he_r2_v248(:,1), he_r3_v248(:,1), he_r4_v248(:,1)];
he_v284.cur = [he_r1_v248(:,2), he_r2_v248(:,2), he_r3_v248(:,2), he_r4_v248(:,2)];
clear he_r1_v248 he_r2_v248 he_r3_v248 he_r4_v248

[he_r1_v248_peaks, r1_loc, r1_width, r1_prom] = findpeaks(abs(he_v284.cur(:,1)));
rem = (1:8);he_r1_v248_peaks(rem) = []; r1_loc(rem) = []; r1_width(rem) = []; r1_prom(rem) = [];

[he_r2_v248_peaks, r2_loc, r2_width, r2_prom] = findpeaks(abs(he_v284.cur(:,2)));
rem = (1:2);he_r2_v248_peaks(rem) = []; r2_loc(rem) = []; r2_width(rem) = []; r2_prom(rem) = [];

[he_r3_v248_peaks, r3_loc, r3_width, r3_prom] = findpeaks(abs(he_v284.cur(:,3)));
rem = (1:2);he_r3_v248_peaks(rem) = []; r3_loc(rem) = []; r3_width(rem) = []; r3_prom(rem) = [];
rem = (3);he_r3_v248_peaks(rem) = []; r3_loc(rem) = []; r3_width(rem) = []; r3_prom(rem) = [];

[he_r4_v248_peaks, r4_loc, r4_width, r4_prom] = findpeaks(abs(he_v284.cur(:,4)));
rem = (1:2);he_r4_v248_peaks(rem) = []; r4_loc(rem) = []; r4_width(rem) = []; r4_prom(rem) = [];
rem = (3);he_r4_v248_peaks(rem) = []; r4_loc(rem) = []; r4_width(rem) = []; r4_prom(rem) = [];

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
he_v284.peaks = [he_r1_v248_peaks, he_r2_v248_peaks, he_r3_v248_peaks,he_r4_v248_peaks];
he_v284.ploc = he_v284.vol([r1_loc, r2_loc, r3_loc, r4_loc]);
he_v284.pwid = [r1_width, r2_width, r3_width, r4_width];
he_v284.pprom = [r1_prom, r2_prom, r3_prom, r4_prom];
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

clear he_r1_v248_peaks he_r2_v248_peaks he_r3_v248_peaks he_r4_v248_peaks r1_loc r2_loc r3_loc r4_loc r1_width r2_width r3_width r4_width r1_prom r2_prom r3_prom r4_prom


he_r1_v196 = load(strcat(path,'data\DAY_2_RUN1_V196.txt'));
he_r2_v196 = load(strcat(path,'data\DAY_2_RUN2_V196.txt'));
he_r3_v196 = load(strcat(path,'data\DAY_2_RUN3_V196.txt'));
he_r4_v196 = load(strcat(path,'data\DAY_2_RUN4_V196.txt'));

he_v196.vol = [he_r1_v196(:,1), he_r2_v196(:,1), he_r3_v196(:,1), he_r4_v196(:,1)];
he_v196.cur = [he_r1_v196(:,2), he_r2_v196(:,2), he_r3_v196(:,2), he_r4_v196(:,2)];
clear he_r1_v196 he_r2_v196 he_r3_v196 he_r4_v196

% he_r1_v196_peaks = findpeaks(abs(he_v196.cur(:,1)));%he_r1_v248_peaks(1:8) = [];
% he_r2_v196_peaks = findpeaks(abs(he_v196.cur(:,2)));%he_r2_v248_peaks(1:2) = [];
% he_r3_v196_peaks = findpeaks(abs(he_v196.cur(:,3)));%he_r3_v248_peaks(1:2) = [];he_r3_v248_peaks(3) = [];
% he_r4_v196_peaks = findpeaks(abs(he_v196.cur(:,4)));%he_r4_v248_peaks(1:2) = [];he_r4_v248_peaks(3) = [];
% he_v196.peaks = [he_r1_v196_peaks, he_r2_v196_peaks, he_r3_v196_peaks,he_r4_v196_peaks];

B = load(strcat(path,'data\DAY_2_Reverse_Bias_V249.txt'));
he_rb.vol = B(:,1); 
he_rb.cur = B(:,2);
clear B
%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

%+++++ Voltage 248 ++++++

he_v284.mean_cur = mean(he_v284.cur,2);
he_v284.std_cur = std(he_v284.cur,0,2); % standard deviation for each point
he_v284.std_cur_all = std(he_v284.mean_cur); %Standard Deviation of the mean.
he_v284.peaks_y = [mean(he_v284.peaks,2), std(he_v284.peaks,0,2)]; 
he_v284.peaks_x = [mean(he_v284.ploc,2), std(he_v284.ploc,0,2)]; 
he_v284.peaks_x(4,1) = he_v284.peaks_x(4,1) - .12;
he_v284.peaks_x(3,1) = he_v284.peaks_x(3,1) - .12;
he_v284.peaks_x(2,1) = he_v284.peaks_x(2,1) - .06;
xneg(1:4) = [0 .6 .12 .12]; xpos = [0 .6 .12 .12];



figure('name', 'Helium Excitation')

subplot(2,1,1)
title('Helium Excitation, Voltage = 1.196V ')
hold on
plot(he_v196.vol(:,1),abs(he_v196.cur(:,1)),'y')
plot(he_v196.vol(:,2),abs(he_v196.cur(:,2)),'r')
plot(he_v196.vol(:,3),abs(he_v196.cur(:,3)),'b')
plot(he_v196.vol(:,4),abs(he_v196.cur(:,4)),'m')
plot(he_v196.vol(:,4),abs(mean(he_v196.cur,2)),'k-.','LineWidth',2)
%errorbar(he_v284.peaks_x(:,1),-he_v284.peaks_y(:,1),-he_v284.peaks_y(:,2),'gx')
hold off
xlabel('Anode Voltage (V)');ylabel('|Detector Current (AU)|')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean')
xlim([19 25.5])


subplot(2,1,2)
title('Helium Excitation, Voltage = 2.84 Volts ')
hold on
plot(he_v284.vol(:,1),abs(he_v284.cur(:,1)),'y')
plot(he_v284.vol(:,2),abs(he_v284.cur(:,2)),'r')
plot(he_v284.vol(:,3),abs(he_v284.cur(:,3)),'b')
plot(he_v284.vol(:,4),abs(he_v284.cur(:,4)),'m')
plot(he_v284.vol(:,4),abs(he_v284.mean_cur),'k-.','LineWidth',2)
%errorbar(he_v284.peaks_x(:,1),-he_v284.peaks_y(:,1),-he_v284.peaks_y(:,2),'gx')
hold off
xlabel('Anode Voltage (V)');ylabel('|Detector Current (AU)|')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean')
xlim([19 25.5])


subplot(2,1,2)
hold on
%plot(he_v284.vol(:,1),he_v284.mean_cur)
%errorbar(he_v284.peaks_x(:,1),-he_v284.peaks_y(:,1),-he_v284.peaks_y(:,2),'rx')
plot(he_rb.vol,he_rb.cur,'k')
hold off
title('Helium Reverse-Bias, Voltage = 2.84 Volts ')
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
xlim([23 30])

figure('Name', 'Helium Excitation - Inverted')
title('Helium Excitation - Inverted')
hold on
plot(he_v284.vol(:,1),abs(he_v284.cur(:,1)),'y')
plot(he_v284.vol(:,2),abs(he_v284.cur(:,2)),'r')
plot(he_v284.vol(:,3),abs(he_v284.cur(:,3)),'b')
plot(he_v284.vol(:,4),abs(he_v284.cur(:,4)),'m')
plot(he_v284.vol(:,4),abs(he_v284.mean_cur),'k-.','LineWidth',2)
plot(he_v284.peaks_x, he_v284.peaks_y,'kv','MarkerSize',9,'MarkerFaceColor',[0.5,0.5,0.5])
%plot(he_rb.vol,he_rb.cur,'k')
%plot(test2(:,1),test2(:,2),'g')
errorbar(he_v284.peaks_x(:,1),he_v284.peaks_y(:,1),he_v284.peaks_y(:,2),'r.')
hold off
xlabel('Anode Voltage (V)');ylabel('| Detector Current (AU) |')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean','Peak Values')
xlim([19 25.5])
clear rem xneg xpos

%+++++ Plot Voltage 196 ++++++
he_v196.mean_cur = mean(he_v196.cur,2);
he_v196.std_cur = std(he_v196.mean_cur,0,2);


figure('name', 'Helium Excitation, V = 1.96 Volts')
subplot(2,1,1)
hold on
plot(he_v196.vol(:,1),he_v196.cur(:,1),'r')
plot(he_v196.vol(:,2),he_v196.cur(:,2),'r')
plot(he_v196.vol(:,3),he_v196.cur(:,3),'r')
plot(he_v196.vol(:,4),he_v196.cur(:,4),'r')
plot(he_v196.vol(:,4),he_v196.mean_cur,'k-.','LineWidth',2)
hold off
title('Helium Excitation, V = 1.96 Volts')
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean')
xlim([19 25.5])

% subplot(2,1,2)
% plot(he_v196.vol(:,4),he_v196.mean_cur,'k-.','LineWidth',2)
% title('Averaged Data')
% xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
% xlim([19 25.5])

%++++++++ Conmparison ++++++++++++++++
figure('Name', 'Helium Excitation - Comparison 2.84V to 1.96V')
hold on
plot(he_v284.vol(:,4),abs(he_v284.mean_cur),'b','LineWidth', 1.5)
plot(he_v196.vol(:,4),abs(he_v196.mean_cur),'r','LineWidth', 1.5)
xlabel('Anode Voltage (V)');ylabel('|Detector Current (AU)|')
legend('Data 2.84V','Data 1.96V');xlim([18 25.5])
title('Comparison between 2.48V and 1.96V')
hold off

%% Mercury Excitation
clc;
me_r1_150deg = load(strcat(path,'data\Mercury_Run_150deg_V6004.txt'));
me_r2_150deg = load(strcat(path,'data\Mercury_Run2_150deg_V6004.txt'));
me_r3_150deg = load(strcat(path,'data\Mercury_Run3_150deg_V6004.txt'));
me_r4_150deg = load(strcat(path,'data\Mercury_Run4_150deg_V6004.txt'));

me_150deg.vol = [me_r1_150deg(:,1), me_r2_150deg(:,1), me_r3_150deg(:,1), me_r4_150deg(:,1)];
me_150deg.cur = [me_r1_150deg(:,2), me_r2_150deg(:,2), me_r3_150deg(:,2), me_r4_150deg(:,2)];
clear me_r1_150deg me_r2_150deg me_r3_150deg me_r4_150deg

me_r1_185deg = load(strcat(path,'data\Mercury_Run_185deg_V6004.txt'));
me_r2_185deg = load(strcat(path,'data\Mercury_Run2_185deg_V6004.txt'));
me_r3_185deg = load(strcat(path,'data\Mercury_Run3_185deg_V6004.txt'));
me_r4_185deg = load(strcat(path,'data\Mercury_Run2_185deg_V6004.txt'));

me_185deg.vol = [me_r1_185deg(:,1), me_r2_185deg(:,1), me_r3_185deg(:,1), me_r4_185deg(:,1)];
me_185deg.cur = [me_r1_185deg(:,2), me_r2_185deg(:,2), me_r3_185deg(:,2), me_r4_185deg(:,2)];
clear me_r1_185deg me_r2_185deg me_r3_185deg me_r4_185deg

%$$$$$$$$$$$$$$$$$$$$ Peaks
[me_150deg.peaks(:,1), me_150deg.loc(:,1), me_150deg.width(:,1), me_150deg.prom(:,1)] = findpeaks(abs(me_150deg.cur(:,1)),'MinPeakProminence',.02);
A = me_150deg.loc(:,1);
me_150deg.loc(:,1) = me_150deg.vol(A,1);
[peak, loc, width, prom] = findpeaks(abs(me_150deg.cur(:,2)),'MinPeakProminence',.02);peak(1) = [];loc(1) = [];width(1) = [];prom(1) = [];
me_150deg.peaks(:,2) = peak; me_150deg.loc(:,2) = me_150deg.vol(loc,2); me_150deg.width(:,2) = width; me_150deg.prom(:,2) = prom;

[peak, loc, width, prom] = findpeaks(abs(me_150deg.cur(:,3)),'MinPeakProminence',.02);peak(1) = [];loc(1) = [];width(1) = [];prom(1) = [];
me_150deg.peaks(:,3) = peak; me_150deg.loc(:,3) = me_150deg.vol(loc,3); me_150deg.width(:,3) = width; me_150deg.prom(:,3) = prom;

[peak, loc, width, prom] = findpeaks(abs(me_150deg.cur(:,4)),'MinPeakProminence',.02);%peak(1) = [];loc(1) = [];width(1) = [];prom(1) = [];
me_150deg.peaks(:,4) = peak; me_150deg.loc(:,4) = me_150deg.vol(loc,4); me_150deg.width(:,4) = width; me_150deg.prom(:,4) = prom;

me_150deg.mean_vol = mean(me_150deg.vol,2);
me_150deg.mean_cur = mean(me_150deg.cur,2);

[peak loc width prom]  = findpeaks(abs(me_150deg.mean_cur),'MinPeakProminence',.02);
me_150deg.mean_peak = -peak; me_150deg.mean_loc = me_150deg.mean_vol(loc); me_150deg.mean_width = width; me_150deg.mean_prom = prom;

%$$$$$$$$$$$$$$$$$$$$$$$$$$

me_150deg.cur_mean = mean(me_150deg.cur,2);
me_150deg.cur_std = std(me_150deg.cur,1,2);

me_185deg.cur_mean = mean(me_185deg.cur,2);
me_185deg.cur_std = std(me_185deg.cur_mean);
me_185deg.mean_cur = mean(me_185deg.cur,2);
me_185deg.std_cut = std(me_185deg.cur_mean);


figure('name', 'Mercury Excitation, V = 6.004 Volts')

subplot(2,1,1)
hold on
plot(me_150deg.vol,me_150deg.cur)
plot(me_150deg.vol(:,1),me_150deg.cur_mean,'k-.','LineWidth',2)
st = sprintf('Mercury Lamp @ ~150%cC',char(176));
title(st)
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean')

subplot(2,1,2)
hold on
plot(me_185deg.vol,me_185deg.cur)
plot(me_185deg.vol(:,1),me_185deg.cur_mean,'k-.','LineWidth',2)
st = sprintf('Mercury Lamp @ ~185%cC',char(176));
title(st)
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean')
hold off


error_pm = mean(me_150deg.loc,2);
error_ps = std(me_150deg.peaks,1,2);
error_ps2 = std(error_pm);
figure('name', 'Mercury Excitation - Inverted, V = 6.004 Volts')

hold on
plot(me_150deg.vol(:,1),abs(me_150deg.cur_mean),'b','LineWidth',2)
plot(me_150deg.vol(:,1),abs(me_150deg.cur_mean)./2,'b--','LineWidth',1)
plot(me_185deg.vol(:,1),abs(me_185deg.cur_mean),'r','LineWidth',2)
plot(me_185deg.vol(:,1),abs(me_185deg.cur_mean)./2,'r--','LineWidth',1)
plot(mean(me_150deg.loc,2), mean(me_150deg.peaks,2),'kv','MarkerSize',9,'MarkerFaceColor',[0.5,0.5,0.5])
%errorbar(mean(me_150deg.loc,2), mean(me_150deg.peaks,2),error_ps)
st = sprintf('Mercury Excitation - Inverted \nComparison between 150%sC and 180%sC',char(176),char(176));
title(st)
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
legend('150, 2.6 Volts','150, 1.6 Volts', '180, 2.6 Volts','180, 1.6 Volts','Mean Caculated Peaks')


%% Data manipulation
He284_Ex_my = mean(he_v284.peaks,2);He284_Ex_sy = std(he_v284.peaks,1,2);
He284_Ex_mx = mean(he_v284.ploc,2);He284_Ex_sx = std(he_v284.ploc,1,2);
He284_Ex_wid = mean(he_v284.pwid,2);He284_Ex_swid = std(he_v284.pwid,1,2);
He284_Ex_prom = mean(he_v284.pprom,2);He284_Ex_sprom = std(he_v284.pprom,1,2);
avg = mean(He284_Ex_mx) - He284_Ex_mx(1);

Perc = He284_Ex_wid./He284_Ex_mx


% [he_v196.peaks,he_v196.loc, he_v196.wid, he_v196.prom]  = findpeaks(abs(mean(he_v196.cur,2)), 'MinPeakProminence',.02)
% he_v196.loc = he_v196.vol(he_v196.loc);
% he_v196.pstd = std(abs(he_v196.cur,2))
% He196_Ex_my = mean(he_v196.peaks,2);He196_Ex_sy = std(he_v196.peaks,1,2);
% He196_Ex_mx = mean(he_v196.ploc,2);He196_Ex_sx = std(he_v196.ploc,1,2);
% He196_Ex_wid = mean(he_v196.pwid,2);He196_Ex_swid = std(he_v196.pwid,1,2);
% He196_Ex_prom = mean(he_v196.pprom,2);He196_Ex_sprom = std(he_v196.pprom,1,2);

[me_185deg.peaks, me_185deg.loc, me_185deg.width, me_185deg.prom] = findpeaks(abs(mean(me_185deg.cur,2)),'MinPeakProminence',.02);
me_185deg.loc = me_185deg.vol(me_185deg.loc);
%%
%close all
%---------Presentation
profile off
figure('Name', 'Helium Excitation - Inverted Presentation')
title('Helium Excitation 2.84V - Inverted')
hold on
plot(he_v284.vol(:,1),abs(he_v284.cur(:,1)),'y')
plot(he_v284.vol(:,2),abs(he_v284.cur(:,2)),'r')
plot(he_v284.vol(:,3),abs(he_v284.cur(:,3)),'b')
plot(he_v284.vol(:,4),abs(he_v284.cur(:,4)),'m')
plot(he_v284.vol(:,4),abs(he_v284.mean_cur),'k-.','LineWidth',2)
plot(he_v284.peaks_x, he_v284.peaks_y,'kv','MarkerSize',9,'MarkerFaceColor',[0.5,0.5,0.5])
kl = mean(he_v284.pprom,2);
for i = 1: length(he_v284.peaks_x(1))
line([he_v284.peaks_x(i,1) he_v284.peaks_x(i,1)],[he_v284.peaks_y(i,1), he_v284.peaks_y(i,1)-kl(i)])
end
%plot(he_rb.vol,he_rb.cur,'k')
%plot(test2(:,1),test2(:,2),'g')
errorbar(he_v284.peaks_x(:,1),he_v284.peaks_y(:,1),he_v284.peaks_y(:,2),'r.')
hold off
xlabel('Anode Voltage (V)');ylabel('| Detector Current (AU) |')
legend('Run 1', 'Run 2', 'Run 3', 'Run 4','Mean','Peak Values')
xlim([19 25.5])
%plot(pres.hevol,pres.hecur);
%%
clc
figure(99)
hold on
plot(he_rb.vol,-he_rb.cur,'b')
plot(he_rb.vol,-he_rb.cur,'or')
title('Helium 2.84V - Reverse Bias')
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')


%%
figure(1)
hold on
plot(mean(he_v196.vol,2),he_v196.mean_cur)
plot(mean(he_v284.vol,2),he_v284.mean_cur)
xlabel('Anode Voltage (V)');ylabel('Detector Current (AU)')
legend('1.96V','2.84V');xlim([15 25])

E = [19.7192 20.5250 22.7960 23.6752 19.7231 20.4169 22.6548 23.0757];
EE = [.2001 .2125 .1354 .2441 .08 .0953 .1287 .0258];
for i = 1: length(E)-4
    
EEE(i) = (E(i)*EE(i) + E(i+4)* EE(i+4))./(EE(i) + EE(i+4));
end
ac = [19.8 20.61 22.90 23.41]
EEE = EEE'
P = abs(EEE - ac)./ ((EEE + ac)./2) * 100