clear, clc;
%% B.1
clear, clc, close all
% Solving the lower, center and upper bands for 1, 1/3, 1/12 octave spectra
% Starting from fc = 1000 Hz, using the relation that fu/fl = 2^(1/n) we
% find all the center frequencies from 20 Hz to 20 kHz. than
% fc = 2^(1/2n)fl -> fl = fc/2^(1/2n).
% fu/fl = 2^(1/n) -> fu = fc*2^(1/2n).
% All for any 1/n octave. Sloving this:

n_oct = [1;3; 12];
fc_vec = zeros(1000,length(n_oct));
freq_range = linspace(20,20000,1000);

lower_to_upper = 2.^(1./n_oct);

for i = 1:length(n_oct)
    fc_upper = [1000];
    while max(fc_upper) < max(freq_range)
        fc_upper(end+1) = fc_upper(end)*lower_to_upper(i);
    end
    fc_lower = [1000];
    while min(fc_lower) > min(freq_range)
        fc_lower(end+1) = fc_lower(end)/lower_to_upper(i);
    end
    fc_vec_temp = [flip(fc_lower) fc_upper(2:end-1)];
    fc_vec(1:length(fc_vec_temp),i) = fc_vec_temp;
end


%% B.2
clear;

f = linspace(10,10^4,1000);

a_1 = 12194;
a_2 = 20.6;
a_3 = 107.7;
a_4 = 737.9;

R_a = (f.^4*a_1^2)./((f.^2 + a_2^2).*sqrt((f.^2 + a_3^2).*(f.^2 + a_4^2)).*(f.^2 + a_1^2));

A = 20*log10(R_a) + 2;

figure;
hold on;
plot(f,A);
grid on;
xlabel('f [Hz]')
ylabel('dB(A)')

A_400_Hz = A(f == 400);
A_200_Hz = A(f == 200);
A_800_Hz = A(f == 800);
scatter(200,A_200_Hz,'filled');
scatter(400,A_400_Hz,'filled');
scatter(800,A_800_Hz,'filled');
set(gca,'xscale','log');

saveFigure(gcf, 'dB(A)', {'png','eps'}, 500);


%% B.3
n = [2, 0.5];
for i = 1:length(n)
    DSPL(i) = 60*log10(n(i));
end




%% C.1
clear, close all;
% Our signal is 10*sin(t). Now we need to compute the Fourier transform
% using FFT, and using the analytical solution

A  = 10;                 % amplitude
fs = [10 100];           % sampling frequency [Hz]
dt = 1./fs;
T_tot = 50;              % sampling time period

% Figure manegment
fig = figure(1);
fig.WindowState = 'maximized';
ax_big = axes('Position',[0.13 0.11 0.775 0.815]);
hold(ax_big,'on')
grid(ax_big,'on')

ax_samll = axes('position',[.65 .175 .25 .25]);
box(ax_samll,'on')
hold(ax_samll,'on')
grid(ax_samll,'on')
axis(ax_samll,'tight')

xlabel(ax_big,'t [sec]')
ylabel(ax_big,'f(x)')
grid (ax_big,'on')
title(ax_big,'Original Signal vs. fft and ifft')


for i = 1:length(dt)
    t = 0:dt(i):T_tot;
    X = A*sin(t);
    F = fft(X);
    F_i = ifft(F);
    % figure(1)
    plot(ax_big,t,X,'DisplayName',['real signal for fs = ' num2str(fs(i)) 'Hz'])
    plot(ax_big,t,F_i,'DisplayName',['transformed signal for fs = ' num2str(fs(i)) 'Hz'],'LineStyle','--')


    indexOfInterest = (t < 1.7) & (t > 1.4); % range of t near perturbation
    plot(ax_samll,t(indexOfInterest),X(indexOfInterest)) % plot on new axes
    plot(ax_samll,t(indexOfInterest),F_i(indexOfInterest)) % plot on new axes
    legend(ax_big,'show','Location','southwest')
    saveFigure(gcf, 'Fourier_Signal', {'png','eps'}, 500);

    N       = length(X);
    df      = fs(i)/N;
    f       = (-N/2:N/2-1)*df;          % frequency axis [-fs/2, fs/2)
    F_shift = fftshift(F)/N;            % center zero freq, normalize

    figure(2);
    stem(f, abs(F_shift));
    grid on;
    xlabel('f [Hz]')
    ylabel('|F(f)|')
    title(['FFT vs analytical, fs = ' num2str(fs(i)) ' Hz'])
    hold on

    % analytical impulses at ±f0 (since sin(t) => omega0 = 1 rad/s)
    f0 = 1/(2*pi);                      % ≈ 0.159 Hz
    scatter(f0,5,'filled','r');
    scatter(-f0,5,'filled','r');
    
    xlim([-5 5])

end

figure(2);
saveFigure(gcf, 'Fourier_FFT_vs_analytical', {'png','eps'}, 500);



%% C.2
clear, clc, close all;
load('DJI_Sound_File\Run35_DJI_white_baseline_B2_CW-calibrated.mat')
% sound(MicData(:,1), Fqs);

p_ref=20e-6;

mic_angle = 105;
mic_idx = find(theta == mic_angle);
N = length(MicData(:,mic_idx));

% Using fft, we'll calculate the single-sided spectral density spectrum
% using:
% p'^2 = int(0 -> inf)[Gpp*df]
% and in the discrete form:
% p'^2 = (1/N^2)*sum(Gpp*Df)

df=Fqs/N;
f=(0:N/2-1)*df; % Taking the first positive side without the nyquist frequency
pos_side_idx = 1:N/2;

p_fluct = MicData(:,mic_idx);  % The mic fluctuations are p'
p_s = fft(p_fluct); % Taking fourier transform on the mic data

p_rms_sqr = abs(p_s(1:N/2)).^2;
p_rms_sqr(2:end-1)=2*(p_rms_sqr(2:end-1)); % Taking the positive side and
% multiplying by 2, exexpt the start and nyquist frequency

Gpp = p_rms_sqr/df;

SPL=10*log10(Gpp*df/(p_ref^2));

figure(1)
semilogx(f,SPL)
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
grid on

saveFigure(gcf, 'SPL_Fourier', {'png','eps'}, 500);

% Now using Welch's method:
Df = 1.235;
L = round(Fqs/Df);
window = hann(L);
noverlap = round(L*0.5);

[Ggg_welch,f_welch] = pwelch(p_fluct,window,noverlap,L,Fqs);
SPL_welch=10*log10(Ggg_welch*Df/(p_ref^2));

figure(2)
semilogx(f_welch,SPL_welch)
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
grid on

saveFigure(gcf, 'SPL_Welch', {'png','eps'}, 500);

omega = RPM/60;
f_welch_norm = f_welch/(B*omega);

figure(3)
semilogx(f_welch_norm,SPL_welch)
xlabel('$\frac{F}{B \cdot \omega}$ [Hz]')
ylabel('SPL [dB]')
grid on

saveFigure(gcf, 'SPL_Welch_Norm', {'png','eps'}, 500);

%% C.3
close all
% Using the favorable Welch's method which produces a cleaner plot.
Df = 1.235;
L = round(Fqs/Df);
window = hann(L);
noverlap = round(L*0.5);

[Ggg_welch,f_welch] = pwelch(MicData,window,noverlap,L,Fqs);


SPL_welch = 10*log10(Ggg_welch*Df/(p_ref^2));

omega = RPM/60;
f_welch_norm = f_welch/(B*omega);

th_vec = deg2rad(theta(:).');


[R, TH] = meshgrid(f_welch_norm, th_vec);
R  = R.';
TH = TH.';

% polar -> Cartesian with the shift so theta = 0 -> (1,0)
Y = R.*cos(TH);
X = R.*sin(TH);

figure;
contourf(X, Y, SPL_welch, 10, 'LineStyle', 'none');
axis equal tight
colormap(hot);
c = colorbar;
c.Label.String = 'SPL [dB]';

% Radial-axis label f/f_b
text(-10, 50, '$f_{norm}$','HorizontalAlignment','center');

saveFigure(gcf, 'SPL_vs_theta', {'png','eps'}, 500);

close all









load handel;
sound(y, Fs);