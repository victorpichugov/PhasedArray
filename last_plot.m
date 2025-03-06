%% Параметры КФАР
K = 9;                                 % Количество элементов
M = 1:9;                               % (1, K)
r = 0.54;                              % Радиус кольца (м)
freqs = (300:100:900) * 1e6;           % Частоты (Гц)
snr = -6;                              % Отношение сигнал/шум (дБ)
az_scan_angles = 0:360;                % (1, N)
L = length(az_scan_angles);
el_scan_angles = 0;

min_res_beamscan = zeros(length(freqs), 1);
min_res_capon = zeros(length(freqs), 1);
min_res_music = zeros(length(freqs), 1);

% Количество экспериментов для усреднения
num_experiments = 100;

for freq_idx = 1:length(freqs)
    freq = freqs(freq_idx);
    k_wave = 2 * pi / freq2wavelen(freq);

    true_angle1 = rand() * 180;

    for delta_angle = 1:180
        true_angle2 = true_angle1 + delta_angle;

        resolution_results_beamscan = zeros(num_experiments, 1);
        resolution_results_capon = zeros(num_experiments, 1);
        resolution_results_music = zeros(num_experiments, 1);

        for exp_idx = 1:num_experiments
            a = exp(-1j * k_wave * r * cosd((360 * (M - 1) / K)' - az_scan_angles) * cosd(el_scan_angles));

            S1 = (sign(randn(1, L)) + 1j * sign(randn(1, L))) / sqrt(2); % Сигнал 1
            S2 = (sign(randn(1, L)) + 1j * sign(randn(1, L))) / sqrt(2); % Сигнал 2

            N = (randn(K, L) + 1j * randn(K, L)) / sqrt(2) * db2mag(-snr); % AWGN

            a1 = exp(-1j * k_wave * r * cosd((360 * (M - 1) / K)' - true_angle1) * cosd(el_scan_angles));
            a2 = exp(-1j * k_wave * r * cosd((360 * (M - 1) / K)' - true_angle2) * cosd(el_scan_angles));

            X = a1 * S1 + a2 * S2 + N;

            Rxx = 1/K * (X * X');

            P_beamscan = real(sum((a' * Rxx) .* a.', 2));
            scanpattern1 = sqrt(abs(P_beamscan));
            scanpattern1 = scanpattern1 / max(scanpattern1);

            theta0 = (true_angle1 + true_angle2) / 2;
            P_theta0_beamscan = scanpattern1(round(theta0) + 1);
            P_theta1_beamscan = scanpattern1(round(true_angle1) + 1);

            if P_theta0_beamscan < P_theta1_beamscan
                resolution_results_beamscan(exp_idx) = 1;
            else
                resolution_results_beamscan(exp_idx) = 0;
            end

            P_capon = 1 ./ real(sum(a' .* (Rxx \ a).', 2));
            scanpattern2 = sqrt(abs(P_capon));
            scanpattern2 = scanpattern2 / max(scanpattern2);

            P_theta0_capon = scanpattern2(round(theta0) + 1);
            P_theta1_capon = scanpattern2(round(true_angle1) + 1);

            if P_theta0_capon < P_theta1_capon
                resolution_results_capon(exp_idx) = 1;
            else
                resolution_results_capon(exp_idx) = 0;
            end

            [eigenvals, eigenvects] = privEig(Rxx);
            noise_eigenvects = eigenvects(:, 3:end);
            D = sum(abs((a' * noise_eigenvects)).^2, 2) + eps(1);
            P_music = 1 ./ D;
            scanpattern3 = sqrt(abs(P_music));
            scanpattern3 = scanpattern3 / max(scanpattern3);

            P_theta0_music = scanpattern3(round(theta0) + 1);
            P_theta1_music = scanpattern3(round(true_angle1) + 1);

            if P_theta0_music < P_theta1_music
                resolution_results_music(exp_idx) = 1;
            else
                resolution_results_music(exp_idx) = 0;
            end
        end

        if mean(resolution_results_beamscan) > 0.9 && min_res_beamscan(freq_idx) == 0
            min_res_beamscan(freq_idx) = delta_angle;
        end
        if mean(resolution_results_capon) > 0.9 && min_res_capon(freq_idx) == 0
            min_res_capon(freq_idx) = delta_angle;
        end
        if mean(resolution_results_music) > 0.9 && min_res_music(freq_idx) == 0
            min_res_music(freq_idx) = delta_angle;
        end

        if min_res_beamscan(freq_idx) > 0 && ...
           min_res_capon(freq_idx) > 0 && ...
           min_res_music(freq_idx) > 0
            break;
        end
    end
end

%% Вывод результатов
fprintf('Минимальный разнос для разрешения сигналов:\n');
fprintf('Частота (МГц) | Beamscan | Capon | MUSIC\n');
for freq_idx = 1:length(freqs)
    fprintf('%d           | %.1f°     | %.1f° | %.1f°\n', ...
        freqs(freq_idx) / 1e6, ...
        min_res_beamscan(freq_idx), ...
        min_res_capon(freq_idx), ...
        min_res_music(freq_idx));
end

%% Визуализация зависимости минимального разноса от частоты
figure;
plot(freqs / 1e6, min_res_beamscan, '-o', 'LineWidth', 1.5); hold on;
plot(freqs / 1e6, min_res_capon, '-s', 'LineWidth', 1.5); hold on;
plot(freqs / 1e6, min_res_music, '-d', 'LineWidth', 1.5); hold on;
xlabel('Частота, МГц');
ylabel('Минимальный разнос, градусы');
grid on;
legend('Beamscan', 'Capon (MVDR)', 'MUSIC');
title('Зависимость минимального разноса от частоты');

function [eigenvals, eigenvects] = privEig(Sx)
    [eigenvects, eigenvalsC] = eig(Sx);
    eigenvals = real(diag(eigenvalsC));
    [eigenvals, indx] = sort(eigenvals, 'descend');
    eigenvects = eigenvects(:, indx);
    eigenvals(eigenvals < 0) = 0;
end