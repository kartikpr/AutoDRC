

fs=48000;
tau_inSamples = 1*fs; %time const
c = exp(-1/tau_inSamples);
[x, fl] = audioread("Weare_uncompressed.mp3");
x = x(:,1);


xFilt_x = x_filter(x).^2; %passing through Butterworth 200Hz LPF
Vms_x = zeros(size(xFilt_x)); 
Vms_1x = 0;
for i = 1:size(Vms_x,1)
    Vms_x(i,1) = c*Vms_1x + (1-c).*xFilt_x(i,1); %Eq from Vickers
    Vms_1x = Vms_x(i,1);
end

N = 100;
B = 1/(N)*ones(N,1); %100 point moving average

%%fvtool
%fvtool(B,1);

Vrms_x = sqrt(filter(B,1,Vms_x));
VdB_x = 20*log10(Vrms_x);
 VdB_x = VdB_x(VdB_x~=0 & isfinite(VdB_x));

meanVdB = mean(VdB_x);
ind_soft = find(VdB_x < prctile(VdB_x,20 ));
ind_loud = find(VdB_x > prctile(VdB_x,85 ));




crossFilt = crossoverFilter( ...
    'NumCrossovers',2, ...
    'CrossoverFrequencies',[300,2500], ...
    'CrossoverSlopes',48, ...
    'SampleRate',fl);
visualize(crossFilt);


[y1,y2,y3] = crossFilt(x);
y = [y1,y2,y3];
for i=1:3
%%Calculating Vms
tau_inSamples = 0.035*fs; %time const
c = exp(-1/tau_inSamples);
xFilt = x_filter(y(:,i)).^2; %passing through Butterworth 200Hz LPF
Vms = zeros(size(xFilt)); 
Vms_1 = 0;
for i = 1:size(Vms,1)
    Vms(i,1) = c*Vms_1 + (1-c).*xFilt(i,1); %Eq from Vickers
    Vms_1 = Vms(i,1);
end

N = 100;
B = 1/(N)*ones(N,1); %100 point moving average

%%fvtool
%fvtool(B,1);

Vrms = sqrt(filter(B,1,Vms));
VdB = 20*log10(Vrms);




%Calculating VdB over time
ind_soft = find(VdB < prctile(VdB,20 ));
ind_loud = find(VdB > prctile(VdB,85 ));


x_lin = linspace(0,length(VdB)+1,length(VdB));
figure(4)
VdBtimesoft = zeros(length(VdB),1);
for i=1:length(ind_soft)
VdBtimesoft(ind_soft(i,1)) = VdB(ind_soft(i,1));

end
VdBtimesoft(VdBtimesoft==0)=-100;
plot(x_lin,VdBtimesoft);
hold on;


VdBtimesloud = zeros(length(VdB),1);
for i=1:length(ind_loud)
VdBtimesloud(ind_loud(i,1)) = VdB(ind_loud(i,1));

end
VdBtimesloud(VdBtimesloud==0)=-100;
plot(x_lin,VdBtimesloud);
hold on;





figure(1);
histogram(VdB,100);
title('input');
hold on;

%%dynamics spread
Vrms_sort=sort(VdB);

for i=1:100
    Ps(i,:) = prctile(VdBtimesloud,i,"all") ;
    Pl(i,:) = prctile(VdBtimesoft,i,"all") ;

end

figure(2);
plot(Ps);
title('input');
hold on;
plot(Pl);
title('input');
hold on;



end

%%Vicker LLML
k = 0.85;
u_sum = 0;
for i=1:length(VdB)
    u(i,1) = k^(-VdB(i));
    u_sum = u_sum + u(i,1);
end

L_uncompressed= 0;
for i=1:length(VdB)
   weight(i,1) = u(i,1)/u_sum;
   L_uncompressed = L_uncompressed + weight(i,1).*VdB(i,1);

end



gaindBL = L_compressed - L_uncompressed;
gainL = 10^(gaindBL/20);
x = x.*gainL;




    L = weight(i) * VdB(i);



%%loudness lufs
[momentary,shortTerm,integrated]=lufs('drift.mp3');
shortTerm_sort=sort(momentary);
figure(3);
plot(shortTerm_sort)
hold on;



[x_Lu, fs] = audioread('drift.mp3');
[loudness, LRA] = integratedLoudness(x_Lu,fs);
fprintf('Loudness before normalization: %.1f LUFS\n',loudness)

target = -11;
gaindB = target - loudness;
gain = 10^(gaindB/20);
x_LuNorm = x_Lu.*gain;
audiowrite('Weare_uncompressed_norm.wav',x_LuNorm,44100);











