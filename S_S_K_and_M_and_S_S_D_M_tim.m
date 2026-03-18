function [S1, S2, S3, n1, n2, n3, l1, l2, l3, q1, q2, q3] = S_S_K_and_M_and_S_S_D_M_tim(n)
% S_S_K_AND_M_AND_S_S_D_M_TIM Сравнительный анализ времени выполнения алгоритмов.
% 
% Эта функция сравнивает три метода решения дискретного уравнения Ляпунова:
% 1. Метод Малышевой.
% 2. Метод S-S(D-M) .
% 3. Метод S-S(K).
%
% Входные параметры:
%   n - размерность исследуемой квадратной матрицы.

    % ------------------------------------------------------------
    % Настройки точности и ограничений
    % ------------------------------------------------------------
    w_max = 1e10;      % Пороговое ограничение нормы
    u_max = 1e10;      % Ограничение числа обусловленности
    ip = 1e-10;        % Требуемая точность вычислений

    % Генерация случайной исходной матрицы A
    A_original = rand(n, n);
    [n, ~] = size(A_original);
    I = eye(n);
    
    % Расчет максимального числа итераций m0 для метода Малышевой
    m0 = round(log2(- (1 + w_max) * log(ip / ((2 + 2 * ip) * sqrt(w_max)))));

    % ================================================================
    % МЕТОД 1: Алгоритм Малышевой 
    % ================================================================
    start_time1 = clock;

    A0 = A_original'; 
    B0 = I;
    p_mat = 2*I;
    d_count = 0;
    
    while (max(abs(p_mat*p_mat - p_mat), [], 'all') > ip * max(abs(p_mat), [], 'all') && d_count <= m0/3)
        d_count = d_count + 1;
        for i1 = 1:3
            S_temp = [-1*B0; A0];
            [Q_mat, ~] = qr(S_temp);
            QQ = Q_mat';
            Q1_sub = QQ(n+1:2*n, 1:n);
            Q2_sub = QQ(n+1:2*n, n+1:2*n);
            A0 = Q1_sub * A0;
            B0 = Q2_sub * B0;
        end

        if (cond(A0 - B0) > u_max)
            disp('Метод Малышевой: Собственное значение слишком близко к единичному кругу.');
            break; 
        end
        p_mat = -((A0 - B0) \ B0);
    end         

    inv_su = inv(B0 + A0);
    H1_final = inv_su * inv_su';
    
    S1 = etime(clock, start_time1);   % Время выполнения Метода 1
    
    % Оценка невязок для Метода 1
    n1 = max(abs(p_mat*p_mat - p_mat), [], 'all'); % Самопроектируемость
    n2 = max(abs(p_mat*A_original' - A_original'* p_mat), [], 'all'); % Коммутируемость
    n3 = max(abs(H1_final - A_original'* H1_final*A_original - p_mat'*p_mat + (I-p_mat)'*(I-p_mat)), [], 'all'); % Уравнение Ляпунова

    % ================================================================
    % МЕТОД 2: S-S(D-M) 
    % ================================================================    
    A_qz = A_original;
    start_time2 = clock;

    [Z_fixed, T_fixed] = schur(A_qz, 'complex'); 
    AA_current = T_fixed; 
    
    % Упорядочивание разложения Шура (собственные значения внутри/снаружи круга)
    select = abs(diag(AA_current)) < 1.0;
    [Z_s, AA_s] = ordschur(Z_fixed, AA_current, select);

    k_dim = sum(select);
    q_dim = n - k_dim;
    
    A11 = AA_s(1:k_dim, 1:k_dim);
    A12 = AA_s(1:k_dim, k_dim+1:n);
    A22 = AA_s(k_dim+1:n, k_dim+1:n);
    
    % Решение уравнения Сильвестера для разделения инвариантных подпространств
    M_sys = sylvester(A11, -A22, -A12);

    I_k = eye(k_dim); I_q = eye(q_dim);
    inv_Q1_transform = [I_k, -M_sys; zeros(q_dim, k_dim), I_q] * Z_s';

    w22_sub_qz = M_sys' * M_sys + I_q;
    inv_BB = I_q / A22; 
    w11_sub_qz = I_k;

    % Итерационный процесс для блока X (внутренняя стабильность)
    aa_iter = A11;
    x0_mat = w11_sub_qz;
    x1_mat = x0_mat + (aa_iter') * x0_mat * aa_iter;
    s_count = 0;
    while (k_dim > 0 && max(abs(x1_mat - x0_mat), [], 'all') > ip&& max(abs(x1_mat), [], 'all') < w_max && s_count < 10)
        s_count = s_count + 1;
        for j_loop = 1:3
            aa_iter = aa_iter * aa_iter;
            x0_mat = x1_mat;
            x1_mat = x0_mat + aa_iter' * x0_mat * aa_iter;
        end
    end

    % Итерационный процесс для блока Y (внешняя стабильность)
    bb_iter = inv_BB;
    y0_mat = bb_iter' * w22_sub_qz * bb_iter;
    y1_mat = y0_mat + bb_iter' * y0_mat * bb_iter;
    s_count = 0;
    while (k_dim < n && max(abs(y1_mat - y0_mat), [], 'all') > ip&& max(abs(y1_mat), [], 'all') < w_max && s_count < 10)
        s_count = s_count + 1;
        for j_loop = 1:3
            bb_iter = bb_iter * bb_iter;
            y0_mat = y1_mat;
            y1_mat = y0_mat + bb_iter' * y0_mat * bb_iter;
        end
    end
    
    % Формирование общего решения H2
    mm_zero = zeros(k_dim, n - k_dim);
    zz_block = [x1_mat, mm_zero; mm_zero', y1_mat];
    H1_final_2 = inv_Q1_transform' * zz_block * inv_Q1_transform;
    
    S2 = etime(clock, start_time2); % Время выполнения Метода 2
    
    % Вычисление проектора для оценки невязок Метода 2
    Z1 = Z_s(:, 1:k_dim); Z2 = Z_s(:, k_dim+1:n);
    pr = Z1 * (Z1' - M_sys * Z2');
    l1 = max(abs(pr*pr - pr), [], 'all');
    l2 = max(abs(pr*A_original - A_original* pr), [], 'all');
    l3 = max(abs(H1_final_2 - A_original'* H1_final_2*A_original - pr'*pr + (I-pr)'*(I-pr)), [], 'all');

    % ================================================================
    % МЕТОД 3: S-S(K)
    % ================================================================
    
    start_time3 = clock;

    [Z_fixed, T_fixed] = schur(A_qz, 'complex'); 
         select = abs(diag(T_fixed)) < 1.0;
    [Z_s, AA_s] = ordschur(Z_fixed, T_fixed, select);
          k_dim = sum(select);
          q_dim = n - k_dim;
    A11 = AA_s(1:k_dim, 1:k_dim);
    A12 = AA_s(1:k_dim, k_dim+1:n);
    A22 = AA_s(k_dim+1:n, k_dim+1:n);
    M_sys = sylvester(A11, -A22, -A12);
      I_k = eye(k_dim);
        I_q = eye(q_dim);
    inv_Q1_transform = [I_k, -M_sys; zeros(q_dim, k_dim),I_q] * Z_s';
    
    % Использование прямого решения для треугольных матриц
    if cond(A22) > u_max
       disp([' Точка r = ', num2str(1), ' находится слишком близко к спектру матрицы'])
     end 
        w22_sub_qz = M_sys' * M_sys + I_q;

            inv_BB = I_q / A22; 
            w11_sub_qz = I_k;
        

        % Итерационный процесс Шура для x
         aa_iter = A11 ;
        x0_mat = w11_sub_qz;
        x1_mat=solve_discrete_triangular(aa_iter,x0_mat);
    
        % Итерационный процесс Шура для y
        bb_iter = inv_BB;
        y0_mat = bb_iter' * w22_sub_qz * bb_iter;
       y1_mat=solve_discrete_triangular( bb_iter,y0_mat);
       
     mm_zero = zeros(k_dim, n - k_dim);
    zz_block = [x1_mat, mm_zero; mm_zero', y1_mat];
    H1_final_3 = inv_Q1_transform' * zz_block * inv_Q1_transform;
    
    S3 = etime(clock, start_time3); % Время выполнения Метода 3
    
    % Оценка невязок для Метода 3
    Z1 = Z_s(:, 1:k_dim); 
    Z2 = Z_s(:, k_dim+1:n);
    pr = Z1 * (Z1' - M_sys * Z2');
    q1 = max(abs(pr*pr - pr), [], 'all');
    q2 = max(abs(pr*A_original - A_original* pr), [], 'all');
    q3 = max(abs(H1_final_3 - A_original'* H1_final_3*A_original - pr'*pr + (I-pr)'*(I-pr)), [], 'all');
end