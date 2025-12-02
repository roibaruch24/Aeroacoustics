clear, clc, close all;
current_path = pwd;
%% B.1
% Solving the lower, center and upper bands for 1, 1/3, 1/12 octave spectra
% Starting from fc = 1000 Hz, using the relation that fu/fl = 2^(1/n) we
% find all the center frequencies from 20 Hz to 20 kHz. than
% fc = 2^(1/2n)fl -> fl = fc/2^(1/2n).
% fu/fl = 2^(1/n) -> fu = fc*2^(1/2n).
% All for any 1/n octave. Sloving this:

n_oct = [1;3; 12];
fc_vec = zeros(1000,length(n_oct));
fl_vec = zeros(1000,length(n_oct));
fu_vec = zeros(1000,length(n_oct));
freq_range = linspace(20,20000,1000);

lower_to_upper = 2.^(1./n_oct);

for i = 1:length(n_oct)
    n = n_oct(i);
    fc_upper = 1000;
    fl_upper = 1000*2^(-1/(2*n));
    fu_upper = 1000*2^(1/(2*n));
    while max(fu_upper) < max(freq_range)
        fc_upper(end+1) = fc_upper(end)*lower_to_upper(i);
        fl_upper(end+1) = fc_upper(end)*2^(-1/(2*n));
        fu_upper(end+1) = fc_upper(end)*2^(1/(2*n));
    end
    fc_lower = 1000;
    fl_lower = 1000*2^(-1/(2*n));
    fu_lower = 1000*2^(1/(2*n));
    while min(fl_lower) > min(freq_range)
        fc_lower(end+1) = fc_lower(end)/lower_to_upper(i);
        fl_lower(end+1) = fc_lower(end)*2^(-1/(2*n));
        fu_lower(end+1) = fc_lower(end)*2^(1/(2*n));
    end
    fc_vec_temp = [flip(fc_lower(1:end-1))  fc_upper(2:end-1)];
    fl_vec_temp = [flip(fl_lower(1:end-1))  fl_upper(2:end-1)];
    fu_vec_temp = [flip(fu_lower(1:end-1))  fu_upper(2:end-1)];

    fc_vec(1:numel(fc_vec_temp),i) = fc_vec_temp(:);
    fl_vec(1:numel(fl_vec_temp),i) = fl_vec_temp(:);
    fu_vec(1:numel(fu_vec_temp),i) = fu_vec_temp(:);
end
idx1  = fc_vec(:,1) > 0;   % valid rows for 1-oct
idx3  = fc_vec(:,2) > 0;   % valid rows for 1/3-oct
idx12 = fc_vec(:,3) > 0;   % valid rows for 1/12-oct

Oct_1 = table( (1:sum(idx1))', fl_vec(idx1,1), fc_vec(idx1,1), fu_vec(idx1,1), ...
    'VariableNames', {'Band','fl_Hz','fc_Hz','fu_Hz'});

Oct_3 = table( (1:sum(idx3))', fl_vec(idx3,2), fc_vec(idx3,2), fu_vec(idx3,2), ...
    'VariableNames', {'Band','fl_Hz','fc_Hz','fu_Hz'});

Oct_12 = table( (1:sum(idx12))', fl_vec(idx12,3), fc_vec(idx12,3), fu_vec(idx12,3), ...
    'VariableNames', {'Band','fl_Hz','fc_Hz','fu_Hz'});

mkdir('Oct_tables')
cd("Oct_tables\")
table2latex(Oct_1, 'Oct_1.tex','1-Octave Bands starting at \(f_c = 1000\) Hz');
table2latex(Oct_3, 'Oct_3.tex','1/3 -Octave Bands starting at \(f_c = 1000\) Hz');
table2latex(Oct_12, 'Oct_12.tex','1/12 -Octave Bands starting at \(f_c = 1000\) Hz');
cd(current_path)

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
text(200, A_200_Hz - 10, sprintf('%.2f dB', A_200_Hz), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center', ...
    'FontSize',10);

text(400, A_400_Hz - 10, sprintf('%.2f dB', A_400_Hz), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center', ...
    'FontSize',10);

text(800, A_800_Hz - 10, sprintf('%.2f dB', A_800_Hz), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center', ...
    'FontSize',10);

set(gca,'xscale','log');
legend('dB(A)','400 Hz','200 Hz','800 Hz','Location','best')

saveFigure(gcf, 'dB(A)', {'png'}, 500);


%% B.3
n = [2, 0.5];
for i = 1:length(n)
    DSPL(i) = 60*log10(n(i));
end




%% C.1
clear, close all;
% Our signal is 10*sin(t). Now we need to compute the Fourier transform
% using FFT, and using the analytical solution


Fs = 100;
k_vec = [2 8 1.5 8.5];     % Number of cycles

for i = 1:length(k_vec)
    k = k_vec(i);
    T = 2*pi*k;              % Sampling period [Sec]
    t = 0:1/Fs:(T - 1/Fs);   % Time vector
    N = length(t);
    Df = Fs/N;
    if mod(N,2)==0
        f_vec = (-N/2 : N/2-1) * Df;     % even N
    else
        f_vec = (-(N-1)/2 : (N-1)/2) * Df; % odd N
    end



    A  = 10;                 % amplitude




    Df = Fs/N;

    X = A*sin(t);
    F = fft(X);
    F_i = ifft(F);

    figure
    hold on;

    xlabel('t [sec]')
    ylabel('f(x)')
    grid ('on')
    title(['Original Signal vs. fft and ifft for ' num2str(k) ' cycles'])

    plot(t,X,'DisplayName',['real signal for Df = ' num2str(Df) 'Hz'])
    plot(t,F_i,'DisplayName',['transformed signal for Df = ' num2str(Df) 'Hz'],'LineStyle','--')


    legend('show','Location','southwest')
    saveFigure(gcf, ['Fourier_Signal_Df_' num2str(Df)] , {'png'}, 500);

    N       = length(X);
    f       = (-N/2:N/2-1)*Df;          % frequency axis [-fs/2, fs/2)
    F_shift = fftshift(F)/N;            % center zero freq, normalize

    figure();
    stem(f, abs(F_shift),'filled');
    grid on;
    xlabel('f [Hz]')
    ylabel('|F(f)|')
    title(['FFT vs analytical, Df = ' num2str(Df) ' Hz for ' num2str(k) ' cycels'])
    hold on

    % analytical impulses at ±f0 (since sin(t) => omega0 = 1 rad/s)
    f0 = 1/(2*pi);                      % ≈ 0.159 Hz
    scatter(f0,5,'filled','r');
    scatter(-f0,5,'filled','r');

    xlim([-1 1])
    legend('FFT','Analytical')
    saveFigure(gcf, ['Fourier_FFT_vs_analytical_Df_' num2str(Df)], {'png'}, 500);

end





%% C.2
clear, clc, close all;
load('DJI_Sound_File\Run35_DJI_white_baseline_B2_CW-calibrated.mat')
% sound(MicData(:,1), Fqs);

p_ref=20e-6;

mic_angle = 105;
mic_idx = find(theta == mic_angle);
N = length(MicData(:,mic_idx));

df=Fqs/N;
f=(0:N/2)*df; % Taking the first positive side without the nyquist frequency

p_fluct = MicData(:,mic_idx);  % The mic fluctuations are p'
p_s = fft(p_fluct); % Taking fourier transform on the mic data

Spp = abs(p_s.^2)/(Fqs*N);
% Spp=(p_rms_sqr_norm); % Taking the positive side and
% multiplying by 2, exexpt the start and nyquist frequency

Gpp = 2*Spp(1:N/2+1);

% Parsvel check FFT
p_rms_time = mean(abs(p_fluct.^2));
p_rms_fft = sum(Gpp)*df;         
ParsevalError = (p_rms_fft - p_rms_time);

fprintf('Parseval Error FFT = %.6e\n', ParsevalError);

SPL=10*log10(Gpp*df/(p_ref^2));

figure(1)
semilogx(f,SPL)
xlabel('Frequency [Hz]')
ylabel('SPL [dB]')
grid on

saveFigure(gcf, 'SPL_Fourier', {'png'}, 500);

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

saveFigure(gcf, 'SPL_Welch', {'png'}, 500);

omega = RPM/60;
f_welch_norm = f_welch/(B*omega);

figure(3)
semilogx(f_welch_norm,SPL_welch)
xlabel('$\frac{F}{B \cdot \omega}$ [Hz]')
ylabel('SPL [dB]')
grid on

% Parsvel check FFT
p_rms_time = mean(abs(p_fluct.^2));
p_rms_fft = sum(Ggg_welch*df);            % p_rms^2 in S domain
ParsevalError = (p_rms_fft - p_rms_time);

fprintf('Parseval Error pwelch = %.6e\n', ParsevalError);


saveFigure(gcf, 'SPL_Welch_Norm', {'png'}, 500);

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

saveFigure(gcf, 'SPL_vs_theta', {'png'}, 500);

close all









load handel;
sound(y, Fs);

%% Help Functions

function saveFigure(figHandle, fileBaseName, formats, resolution)
% saveFigure  Save a MATLAB figure in multiple formats
%
%   saveFigure(figHandle, fileBaseName, formats, resolution)
%
%   Inputs:
%     figHandle    – Handle to the figure (e.g. gcf)
%     fileBaseName – Base name of the file (no extension)
%     formats      – Cell array of strings, e.g. {'png','jpg','pdf'}
%     resolution   – DPI resolution (optional, default = 300)
%
%   Example:
%     % Save current figure as PNG, JPEG, and PDF at 300 dpi
%     saveFigure(gcf, 'airfoil_lift_curve', {'png','jpeg','pdf'}, 300);

    if nargin < 4
        resolution = 300;
    end

    % Create "plots" folder in current working directory if it doesn't exist
    plotDir = fullfile(pwd, 'plots');
    if ~exist(plotDir, 'dir')
        mkdir(plotDir);
    end

    for k = 1:numel(formats)
        fmt = lower(formats{k});

        switch fmt
            case 'png'
                fileName = fullfile(plotDir, [fileBaseName '.png']);
                print(figHandle, fileName, '-dpng', ['-r' num2str(resolution)]);

            case {'jpg','jpeg'}
                fileName = fullfile(plotDir, [fileBaseName '.jpg']);
                print(figHandle, fileName, '-djpeg', ['-r' num2str(resolution)]);

            case {'tif','tiff'}
                fileName = fullfile(plotDir, [fileBaseName '.tif']);
                print(figHandle, fileName, '-dtiff', ['-r' num2str(resolution)]);

            case 'bmp'
                fileName = fullfile(plotDir, [fileBaseName '.bmp']);
                print(figHandle, fileName, '-dbmp', ['-r' num2str(resolution)]);

            case 'eps'
                fileName = fullfile(plotDir, [fileBaseName '.eps']);
                print(figHandle, fileName, '-depsc', ['-r' num2str(resolution)]);

            case 'pdf'
                fileName = fullfile(plotDir, [fileBaseName '.pdf']);
                print(figHandle, fileName, '-dpdf', ['-r' num2str(resolution)]);

            otherwise
                warning('saveFigure:unknownFormat', ...
                        'Format "%s" not supported – skipping.', fmt);
        end
    end
end

function table2latex(tbl, filename, caption, label)

    % Optional caption/label
    if nargin < 3, caption = ''; end
    if nargin < 4, label = ''; end

    % Replace _ with \_
    escapeUnderscore = @(s) strrep(s, '_', '\_');

    fid = fopen(filename,'w');

    fprintf(fid, '\\begin{table}[h]\n');
    fprintf(fid, '\\centering\n');

    fprintf(fid, '\\begin{tabular}{%s}\n', repmat('c',1,width(tbl)));
    fprintf(fid, '\\hline\n');

    % Header row
    for i = 1:width(tbl)
        header = escapeUnderscore(tbl.Properties.VariableNames{i});
        fprintf(fid, '\\(%s\\)', header);
        if i < width(tbl), fprintf(fid, ' & '); end
    end
    fprintf(fid, ' \\\\ \\hline\n');

    % Body rows
    for r = 1:height(tbl)
        for c = 1:width(tbl)
            fprintf(fid, '\\(%g\\)', tbl{r,c});
            if c < width(tbl), fprintf(fid, ' & '); end
        end
        fprintf(fid, ' \\\\ \n');
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');

    % Bottom caption + label
    if ~isempty(caption)
        fprintf(fid, '\\caption{%s}\n', caption);
    end
    if ~isempty(label)
        fprintf(fid, '\\label{%s}\n', label);
    end

    fprintf(fid, '\\end{table}\n');

    fclose(fid);
end




