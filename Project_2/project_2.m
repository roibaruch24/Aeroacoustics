clc, clear, close all;
%% A.3
k_a = logspace(-1,1,1000);
NC = 2*exp(-1i*k_a)./(2 - 2i*k_a - k_a.^2);
NC_db = mag2db(abs(NC));

ka_1_idx  = find(abs(NC_db) < 1,    1, 'last');
ka_09_idx = find(abs(NC_db) < 0.09, 1, 'last');

ka_1  = k_a(ka_1_idx);
ka_09 = k_a(ka_09_idx);

kr = logspace(-1,2,1000);
NF = 1 - 1./(1i*kr);
NF_db = mag2db(abs(NF));

kr_1_idx  = find(abs(NF_db) < 1);
kr_09_idx = find(abs(NF_db) < 0.09);

kr_1  = kr(kr_1_idx(1));
kr_09 = kr(kr_09_idx(1));


figure
plot(k_a, abs(NC_db), 'LineWidth', 1.3)
hold on
grid on
set(gca,'XScale','log')

scatter(ka_1,  abs(NC_db(ka_1_idx)),  60, 'filled')
scatter(ka_09, abs(NC_db(ka_09_idx)), 60, 'filled')

text(ka_1, abs(NC_db(ka_1_idx)), ...
    sprintf('1 dB at $ka=%.2f$', ka_1), ...
    'Interpreter','latex', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center')

text(ka_09*0.8, abs(NC_db(ka_09_idx))*1.3, ...
    sprintf('0.09 dB at $ka=%.2f$', ka_09), ...
    'Interpreter','latex')

ylim([0 8])
xlabel('$ka$','Interpreter','latex')
ylabel('dB','Interpreter','latex')
title('Compact Approximation','Interpreter','latex')

saveFigure(gcf, 'Compact_Approximation', {'png'}, 500);


figure
plot(kr, abs(NF_db), 'LineWidth', 1.3)
hold on
grid on
set(gca,'XScale','log')

scatter(kr_1,  abs(NF_db(kr_1_idx(1))),  60, 'filled')
scatter(kr_09, abs(NF_db(kr_09_idx(1))), 60, 'filled')

text(kr_1, abs(NF_db(kr_1_idx(1))), ...
    sprintf('1 dB at $kr=%.2f$', kr_1), ...
    'Interpreter','latex', ...
    'VerticalAlignment','bottom')

text(kr_09*1.2, abs(NF_db(kr_09_idx(1)))*1.3, ...
    sprintf('0.09 dB at $kr=%.2f$', kr_09), ...
    'Interpreter','latex')

ylim([0 8])
xlabel('$kr$','Interpreter','latex')
ylabel('dB','Interpreter','latex')
title('Far Field','Interpreter','latex')

saveFigure(gcf, 'Far_Field', {'png'}, 500);

%% B.2.1
M_vec = [0.2, 0.4, 0.6];
ctau_d = linspace(-20, 20, 100);

ct_d       = zeros(length(M_vec), length(ctau_d));
omega_norm = zeros(length(M_vec), length(ctau_d));

for i = 1:length(M_vec)
    ct_d(i,:) = ctau_d + sqrt(M_vec(i)^2 .* ctau_d.^2 + 1);
    omega_norm(i,:) = 1 ./ (1 + M_vec(i)^2 .* ctau_d ./ sqrt(M_vec(i)^2 .* ctau_d.^2 + 1));
end

figure;hold on;grid on;
for i = 1:length(M_vec)
    plot(ct_d(i,:), omega_norm(i,:),'DisplayName',['M = ' num2str(M_vec(i))])
end

xlim([-12 12])
xlabel('$\frac{c\tau_d}{d}$')
ylabel('$\frac{\omega}{\omega_0}$')

legend Location best

saveFigure(gcf, 'B_2_1', {'png'}, 500);
%% C.4.3

k = 0.57;
M_vec = 0.1:0.1:0.8;
Rossiter_n_modes = 1:6;

LoD = [4 6 8 10];
alpha_v = [0.25 0.38 0.54 0.58];


figure; grid on; hold on;
for n = Rossiter_n_modes
    St = (n - alpha_v(1)) ./ (1/k + M_vec);
    plot(M_vec, St, 'DisplayName', ['n = ' num2str(n)]);
end

xlabel('$M_\infty$','Interpreter','latex');
ylabel('$St_L$','Interpreter','latex');
title(['Rossiter modes, $LoD = ' num2str(LoD(1)) '$'], ...
    'Interpreter','latex');
legend('Location','best');
saveFigure(gcf, 'StL_M', {'png'}, 500);

figure; grid on; hold on;
for n = Rossiter_n_modes
    St = (n - alpha_v) ./ (1/k + M_vec(3));
    plot(LoD, St,'DisplayName', ['n = ' num2str(n)]);
end

xlabel('$LoD$','Interpreter','latex');
ylabel('$St_L$','Interpreter','latex');
title(['Rossiter modes vs $LoD$, $M_\infty = ' num2str(M_vec(3)) '$'], ...
      'Interpreter','latex');
legend('Location','best');
saveFigure(gcf, 'StL_Lod', {'png'}, 500);

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

