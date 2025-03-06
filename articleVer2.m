clc; clear; close all;

%% Параметры КФАР
K = 9;                                 % Количество элементов
M = 1:9;                               % (1, K)
r = 0.54;                              % Радиус кольца (м)
% freqs = [300, 600, 900]*1e6; % Частоты (Гц)
freqs = 900 *1e6;                      % Частота (Гц)
k_wave = 2 * pi / freq2wavelen(freqs); % волновое число

%% Источник сигнала
snr = 3;                   % Отношение сигнал/шум (дБ)
true_angle = 50;           % Истинный угол (градусы)
az_scan_angles = 0:360;    % (1, N)
L = length(az_scan_angles);
el_scan_angles = 0;

% Вектор направленности (Steering Vectors) (K, L)
a = exp(-1j * k_wave * r * cosd((360 * (M - 1) / K)' - az_scan_angles) * cosd(el_scan_angles));

% Генерация принятого сигнала
% S = (sign(randn(K, L)) + 1j*sign(randn(K, L))) / sqrt(2); % QPSK     % Сигнал источника QPSK
S = (sign(randn(1, L)) + 1j*sign(randn(1, L))) / sqrt(2); % QPSK     % Сигнал источника QPSK
N = (randn(K, L) + 1j*randn(K, L)) / sqrt(2) * db2mag(-snr); % AWGN
% Направляющий вектор с направления откуда идет сигнал
a1=exp(-1j * k_wave * r * cosd((360 * (M - 1) / K)' - true_angle) * cosd(el_scan_angles));

% X = a .* S + N; % (K, L)
X=a1*S+N;

% Оценка корреляционной матрицы 1/L
Rxx = 1/K*(X * X');       % (K, K)

%% Метод сканирования (Beamscan)
P_beamscan = real(sum((a'*Rxx).*a.',2)); % (L, 1) 
scanpattern1 = sqrt(abs(P_beamscan));     % (L, 1)
scanpattern1 = scanpattern1 / max(scanpattern1);


%% Метод Capon (MVDR)
P_capon = 1./real(sum(a'.*(Rxx\a).',2)); % (L, 1) 
scanpattern2 = sqrt(abs(P_capon));     % (L, 1)
scanpattern2 = scanpattern2 / max(scanpattern2);

%% Метод MUSIC (Subspace)

[eigenvals, eigenvects] = privEig(Rxx);
noise_eigenvects = eigenvects(:, 2:end);
D = sum(abs((a'*noise_eigenvects)).^2,2)+eps(1);
P_music = 1./D;
scanpattern3 = sqrt(abs(P_music));     % (L, 1)
scanpattern3 = scanpattern3 / max(scanpattern3);


% Расчет для Beamscan
main_lobe_width1 = compute_mainlobe_width(scanpattern1);
sidelobe_level1 = compute_sidelobe_level(scanpattern1);
% Расчет для Capon
main_lobe_width2 = compute_mainlobe_width(scanpattern2);
sidelobe_level2 = compute_sidelobe_level(scanpattern2);

% Расчет для MUSIC
main_lobe_width3 = compute_mainlobe_width(scanpattern3);
sidelobe_level3 = compute_sidelobe_level(scanpattern3);

% Вывод результатов
fprintf('Рабочая частота : %.1f МГц\n', freqs / 1e6);
fprintf('Beamscan: Ширина ГЛ = %.1f°, УБЛ = %.1f%%\n', main_lobe_width1, sidelobe_level1);
fprintf('Capon: Ширина ГЛ = %.1f°, УБЛ = %.1f%%\n', main_lobe_width2, sidelobe_level2);
fprintf('MUSIC: Ширина ГЛ = %.1f°, УБЛ = %.1f%%\n', main_lobe_width3, sidelobe_level3);
fprintf('------------------------------------\n');

%% Визуализация
figure;
plot(az_scan_angles, scanpattern1, 'r', 'LineWidth', 1.5); hold on;
plot(az_scan_angles, scanpattern2, 'g', 'LineWidth', 1.5); hold on;
plot(az_scan_angles, scanpattern3, 'b', 'LineWidth', 1.5); hold on;
xlabel('Угол (градусы)'); 
grid on;
legend('Beamscan', 'Capon (MVDR)', 'MUSIC');
title(sprintf('Метод сканирования на частоте %d МГц', freqs/1e6));

%% Диаграммы направленности в полярных координатах
figure;
subplot(1, 3, 1);
% polarplot(scanpattern1, 'r', 'LineWidth', 1.5); hold on;
polarpattern(scanpattern1);
title('Beamscan');

subplot(1, 3, 2);
% polarplot(scanpattern2, 'g', 'LineWidth', 1.5); hold on;
polarpattern(scanpattern2);
title('Capon (MVDR)');

subplot(1, 3, 3);
% polarplot(scanpattern3, 'b', 'LineWidth', 1.5); hold on;
polarpattern(scanpattern3);
title('MUSIC');


%%  Численные характеристики
% Обновление переменных
freqs = (300:100:900) * 1e6; % Частоты (Гц)
snr = -6;                              % Отношение сигнал/шум (дБ)
num_experiments = 100;                 % Количество экспериментов
k_wave = 2 * pi ./ freq2wavelen(freqs); % Волновое число
lambda = freq2wavelen(freqs);
d = 3e8 ./ (2 * freqs);

mse_table = zeros(length(freqs), 3);
errors = zeros(num_experiments, 3); % Ошибки

for idx = 1:length(freqs)   
    true_angle = rand * 360;
  
    for exp_idx = 1:num_experiments
       
        a = exp(-1j * k_wave(idx) * r * cosd((360 * (M - 1) / K)' - az_scan_angles) * cosd(el_scan_angles));
        
        S = (sign(randn(1, L)) + 1j * sign(randn(1, L))) / sqrt(2);
        N = (randn(K, L) + 1j * randn(K, L)) / sqrt(2) * db2mag(-snr);
        a1 = exp(-1j * k_wave(idx) * r * cosd((360 * (M - 1) / K)' - true_angle) * cosd(el_scan_angles));
        X = a1 * S + N;
        
        Rxx = 1/K * (X * X');
        
        % Метод Beamscan
        P_beamscan = real(sum((a' * Rxx) .* a.', 2));
        scanpattern1 = sqrt(abs(P_beamscan));
        scanpattern1 = scanpattern1 / max(scanpattern1);
        [~, est_angle1] = max(scanpattern1);
        %
        if abs(az_scan_angles(est_angle1) - true_angle) > 3
            scanpattern1(est_angle1) = 0;
            [~, est_angle1] = max(scanpattern1);
        end
        % 
        errors(exp_idx, 1) = (az_scan_angles(est_angle1) - true_angle)^2;
        
        % Метод Capon
        P_capon = 1 ./ real(sum(a' .* (Rxx \ a).', 2));
        scanpattern2 = sqrt(abs(P_capon));
        scanpattern2 = scanpattern2 / max(scanpattern2);
        [~, est_angle2] = max(scanpattern2);
        %
        if abs(az_scan_angles(est_angle2) - true_angle) > 3
            scanpattern2(est_angle2) = 0;
            [~, est_angle2] = max(scanpattern2);
        end
        %
        errors(exp_idx, 2) = (az_scan_angles(est_angle2) - true_angle)^2;
        
        % Метод MUSIC
        [eigenvals, eigenvects] = privEig(Rxx);
        noise_eigenvects = eigenvects(:, 3:end);
        D = sum(abs((a' * noise_eigenvects)).^2, 2) + eps(1);
        P_music = 1 ./ D;
        scanpattern3 = sqrt(abs(P_music));
        scanpattern3 = scanpattern3 / max(scanpattern3);
        [~, est_angle3] = max(scanpattern3);
        %
        if abs(az_scan_angles(est_angle3) - true_angle) > 3
            scanpattern3(est_angle3) = 0;
            [~, est_angle3] = max(scanpattern3);
        end
        %
        errors(exp_idx, 3) = (az_scan_angles(est_angle3) - true_angle)^2;
    end

    % Расчет СКО для текущей частоты
    mse_table(idx, 1) = sqrt(mean(errors(:, 1)));
    mse_table(idx, 2) = sqrt(mean(errors(:, 2)));
    mse_table(idx, 3) = sqrt(mean(errors(:, 3)));
end

%%  Зависимость СКО пеленга от рабочей частоты
fprintf('СКО пеленга для SNR = -6 дБ:\n');
fprintf('Частота, МГц | Beamscan | MVDR | MUSIC\n');
for idx = 1:length(freqs)
    fprintf('%d          | %.4f°          | %.4f°       | %.4f°\n', ...
        freqs(idx)/1e6, mse_table(idx, 1), mse_table(idx, 2), mse_table(idx, 3));
end

figure;
plot(freqs./1e6, mse_table(:, 1), ...
    freqs./1e6, mse_table(:, 2), '--', ...
    freqs./1e6, mse_table(:, 3), ':', LineWidth=1.5);
xlabel('Частота, МГц'); 
ylabel('СКО пеленга, градусы'); 
grid on;
legend('Beamscan', 'Capon (MVDR)', 'MUSIC');
title(sprintf('Зависимость СКО пеленга от рабочей частоты'));
fprintf('------------------------------------\n');


%%  Зависимость минимального разноса между углами от частоты
fprintf('Разрешающая способность алгоритмов\n');
fprintf('Частота, МГц | Beamscan | MVDR | MUSIC\n');
for idx = 1:length(freqs)
    fprintf('%d          | %.4f°          | %.4f°       | %.4f°\n', ...
        freqs(idx)/1e6, res_table(idx, 1), res_table(idx, 2), res_table(idx, 3));
end

figure;
plot(freqs./1e6, res_table(:, 1), ...
    freqs./1e6, res_table(:, 2), '--', ...
    freqs./1e6, res_table(:, 3), ':', LineWidth=1.5);
xlabel('Частота, МГц'); 
ylabel('Разнос между углами, градусы'); 
grid on;
legend('Beamscan', 'Capon (MVDR)', 'MUSIC');
title(sprintf(' Зависимость минимального разноса между углами от частоты'));
fprintf('------------------------------------\n');


%% Функции
function [eigenvals, eigenvects] = privEig(Sx)
    [eigenvects, eigenvalsC] = eig(Sx);
    eigenvals = real(eigenvalsC);
    [eigenvals,indx] = sort(diag(eigenvals),'descend');
    eigenvects= eigenvects(:,indx);
    eigenvals(eigenvals<0) = 0;
end

%% Функция для проверки разрешения
function mainlobe_width = compute_mainlobe_width(spectrum)
    abs_spectrum = abs(spectrum);
    [~, max_idx] = max(abs_spectrum);
    
    half_power = abs_spectrum(max_idx) / sqrt(2);
    
    left_idx = max_idx;
    while left_idx > 1 && abs_spectrum(left_idx) > half_power
        left_idx = left_idx - 1;
    end
    
    right_idx = max_idx;
    while right_idx < length(abs_spectrum) && abs_spectrum(right_idx) > half_power
        right_idx = right_idx + 1;
    end
    
    mainlobe_width = right_idx - left_idx;
end

function sidelobe_level = compute_sidelobe_level(spectrum)
    abs_spectrum = abs(spectrum);
    [~, max_idx] = max(abs_spectrum);
    
    mainlobe_width = compute_mainlobe_width(spectrum);

    left_bound = max(1, floor(max_idx - mainlobe_width / 2));
    right_bound = min(length(abs_spectrum), ceil(max_idx + mainlobe_width / 2));


    sidelobe_region = [1:left_bound-1, right_bound+1:length(abs_spectrum)];

    sidelobe_max = max(abs_spectrum(sidelobe_region));
    
    sidelobe_level = (sidelobe_max / max(abs_spectrum)) * 100;
end