function [] = tim_c_malshv(n)
% TIM_C_MALSHV Основной скрипт для проведения сравнительного анализа методов.
%
% Функция выполняет серию тестов для матриц различной размерности,
% собирает статистику по времени выполнения и точности (невязкам),
% сохраняет результаты и строит графики.
%
% Входные параметры:
%   n - Максимальная размерность матрицы для тестирования.

    % --- Настройки эксперимента ---
    k = 10:10:n;       % Вектор размерностей матриц (шаг 10)
    c = 15;            % Количество повторений для каждой размерности (для усреднения)
    m = length(k);     % Общее количество шагов по размерности
    
    % Инициализация ячеек для хранения результатов
    z1 = cell(1,m); z2 = cell(1,m); z22 = cell(1,m);
    z3 = cell(1,m); z4 = cell(1,m); z5 = cell(1,m);
    z6 = cell(1,m); z7 = cell(1,m); z8 = cell(1,m);
    z9 = cell(1,m); z10 = cell(1,m); z11 = cell(1,m);

    % --- Цикл проведения вычислений ---
    for j = 1:m
        fprintf('Тестирование размерности n = %d...\n', k(j));
        
        % Временные массивы для итераций внутри одной размерности
        ww1 = zeros(1,c); ww2 = zeros(1,c); ww3 = zeros(1,c);
        uu1 = zeros(1,c); uu2 = zeros(1,c); uu3 = zeros(1,c);
        vv1 = zeros(1,c); vv2 = zeros(1,c); vv3 = zeros(1,c);
        gg1 = zeros(1,c); gg2 = zeros(1,c); gg3 = zeros(1,c);

        for j1 = 1:c
            % Вызов основной функции сравнения трех методов
            [w1, w2, w3, u1, u2, u3, v1, v2, v3, g1, g2, g3] = S_S_K_and_M_and_S_S_D_M_tim(k(j));
            
            % Сбор данных о времени
            ww1(j1) = w1; ww2(j1) = w2; ww3(j1) = w3;
            % Сбор данных о невязках
            uu1(j1) = u1; uu2(j1) = u2; uu3(j1) = u3;
            vv1(j1) = v1; vv2(j1) = v2; vv3(j1) = v3;
            gg1(j1) = g1; gg2(j1) = g2; gg3(j1) = g3;
        end
        
        % Сохранение итерационных данных в cell-массивы
        z1{j}=ww1; z2{j}=ww2; z22{j}=ww3;
        z3{j}=uu1; z4{j}=uu2; z5{j}=uu3;
        z6{j}=vv1; z7{j}=vv2; z8{j}=vv3;
        z9{j}=gg1; z10{j}=gg2; z11{j}=gg3;
    end

    % --- Сохранение результатов в файл ---
    save_filename = 'Simulation_Results.mat';
    save(save_filename, 'k', 'c', 'z1', 'z2', 'z22', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10', 'z11');
    fprintf('Результаты успешно сохранены в файл: %s\n', save_filename);

    % ================================================================
    % ПОСТРОЕНИЕ ГРАФИКОВ
    % ================================================================
    
    % --- График 1: Сравнение времени выполнения (все методы) ---
    figure; hold on;
    mean_z1 = cellfun(@mean, z1);
    mean_z2 = cellfun(@mean, z2);
    mean_z22 = cellfun(@mean, z22);

    for j = 1:length(k)
        plot(ones(1,c)*k(j), z1{j}, 'b+', 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), z2{j}, 'ro', 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), z22{j}, 'k*', 'MarkerSize', 4, 'HandleVisibility', 'off');
    end

    plot(k, mean_z1, 'b-', 'LineWidth', 2);
    plot(k, mean_z2, 'r-', 'LineWidth', 2);
    plot(k, mean_z22, 'k-', 'LineWidth', 2);

    xlabel('Размерность матрицы, n'); ylabel('Время (сек)');
    legend('Метод М','Метод S-S(D-M)','Метод S-S(K)', 'Location', 'northwest');
     set(gca, 'FontSize', 25);

    
        figure; hold on;
      for j = 1:length(k)
        plot(ones(1,c)*k(j), z1{j}, 'b+', 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), z2{j}, 'ro', 'MarkerSize', 5, 'HandleVisibility', 'off');
    end

    plot(k, mean_z1, 'b-', 'LineWidth', 2);
    plot(k, mean_z2, 'r-', 'LineWidth', 2);

    xlabel('Размерность матрицы, n'); ylabel('Время (сек)');
    legend('Метод М','Метод S-S(D-M)', 'Location', 'northwest');
     set(gca, 'FontSize', 25);

    % --- График 2: Точность выполнения условия проектора ||P^2 - P|| ---
    figure; hold on;
    mean_z3 = cellfun(@mean, z3);
    mean_z6 = cellfun(@mean, z6);
    mean_z9 = cellfun(@mean, z9);

    for j = 1:length(k)
        plot(ones(1,c)*k(j), log10(z3{j}), 'b+', 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), log10(z6{j}), 'ro', 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), log10(z9{j}), 'k*', 'HandleVisibility', 'off');
    end
    plot(k, log10(mean_z3), 'b-', 'LineWidth', 2);
    plot(k, log10(mean_z6), 'r-', 'LineWidth', 2);
    plot(k, log10(mean_z9), 'k-', 'LineWidth', 2);

    xlabel('Размерность матрицы, n'); ylabel('log10 ||P^2 - P||');
    legend('Метод М','Метод S-S(D-M)','Метод S-S(K)', 'Location', 'northwest');
     set(gca, 'FontSize', 25);

    % --- График 3: Невязка матричного уравнения ---
    figure; hold on;
    mean_z5 = cellfun(@mean, z5);
    mean_z8 = cellfun(@mean, z8);
    mean_z11 = cellfun(@mean, z11);

    for j = 1:length(k)
        plot(ones(1,c)*k(j), log10(z5{j}), 'b+', 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), log10(z8{j}), 'ro', 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), log10(z11{j}), 'k*', 'HandleVisibility', 'off');
    end
    plot(k, log10(mean_z5), 'b-', 'LineWidth', 2);
    plot(k, log10(mean_z8), 'r-', 'LineWidth', 2);
    plot(k, log10(mean_z11), 'k-', 'LineWidth', 2);

    xlabel('Размерность матрицы, n'); ylabel('Невязка матричного уравнения');
    legend('Метод М','Метод S-S(D-M)','Метод S-S(K)', 'Location', 'northwest');
     set(gca, 'FontSize', 25);



    figure;
    hold on;
    mean_z4 = cellfun(@mean, z4);
    mean_z7 = cellfun(@mean, z7);
    mean_z10 = cellfun(@mean, z10);

    for j = 1:length(k)
        plot(ones(1,c)*k(j),log10( z4{j}), 'b+', 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), log10( z7{j}), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), log10( z10{j}), 'k*', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    end
     plot(k, log10(mean_z4), 'b-', 'LineWidth', 2, 'MarkerSize', 10);
     plot(k, log10(mean_z7), 'r-', 'LineWidth', 2, 'MarkerSize', 9, 'MarkerFaceColor', 'r');
     plot(k, log10(mean_z10), 'k-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    xlabel('Размерность матрицы, n');
    ylabel('$\|P*A-A*P\|$', 'Interpreter', 'latex');

     legend('Метод М','Метод S-S(D-M)','Метод S-S(K)', 'Location', 'northwest');

    xticks(k);
    set(gca, 'FontSize', 25);