%% Практикум по ОУ

% Лабораторная работа №1
% Линейная задача быстродействия
% Вариант 7

%% Входные параметры
clc;
clear;

% Пример 1: вещественные собственные значения
A = [[-2, 0]; [0, 1]]; 
B = [[3, 0]; [0, 4]];
x1 = [4; 2];

% Пример 2: комплексные собственные значения
% A = [[1, -2]; [1, 1]];
% B = [[3, 0]; [0, 4]];
% x1 = [-2; 8];

% Пример 3: вырожденное управление
% A = [[-1, 3]; [-3, 1]]; 
% B = [[0, 1]; [0, 0]];
% x1 = [-2; 2];


% Пример 4: отсутствие непрерывности переменной T по целевому множеству {x_1}
% A = [[3, -15]; [3, 3]]; 
% B = [[0, 1]; [0, 0]];
% x1 = [-1; 1.92];
% x1 = [-1; 1.8];
% a = 2;
% b = 7.2;

a = 2;
b = 16;

f = [0; 0];

t0 = 0;

r = 4;
alpha = 2;

%% Проверка корректности входных параметров

if t0 < 0 || r < 0 || a < 0 || b < 0
    error('Input parameters: t0, r, a, b - are incorrect');
end

%% Проверка условий существования множества X0

 if b < a^2
    error('Input parameters: X0 is an empty set because b < a^2');
 end
 
%% Проверка вложенности X1 в X0

is_in_X0 = @ (x, y) ( ...
        (x >= 0 && x <= sqrt(b) - a && abs(y) <= sqrt(b - (x + a)^2)) ...
        || ... 
        (x >= a - sqrt(b) && x <= 0 && abs(y) <= sqrt(b - (x - a)^2)) ...
        );
    
    
 if is_in_X0(x1(1), x1(2)) == true
     disp('x_1 belongs to X_0');
     disp('Решение линейной задачи быстродействия:');
     disp(['Время быстродействия T = ', num2str(0)]);
     disp(['Точка прибытия в X_0: (', num2str(x1(1)), ', ', num2str(x1(2)), ')' ]);
     disp(['Погрешность выполнения условия трансверсальности в точке прибытия: ', num2str(0)]);
 else

%% Параметры для решения задачи в обратном времени

alpha_n_point_cnt = 12;
alpha_start = 0;
alpha_end = 2 * pi;
alpha_seg_lim = 1e-2;
alpha_opt = pi;
psi_0_idx_opt = 0;

eps = 1e-6;
ode45_max_dt_step = 0.01;

C = -A;
D = -B;
g = -f;

T = 1;
t_opt = -T - 1;

t_grid_point_cnt = 200;
t_grid = linspace(-T, -t0, t_grid_point_cnt);

psi_0 = zeros(2, alpha_n_point_cnt);
psi = zeros(size(psi_0, 2), size(t_grid, 2), 2);
opt_control = zeros(size(psi));
y = zeros(size(psi));

is_solved = true;

scalar_prod = @ (x, y) x(1) * y(1) + x(2) * y(2);

save oc_lw1_temp.mat a b r alpha eps C D g is_in_X0;

%% Решение задачи

% left_circle_x = linspace(a - sqrt(b), 0, 500);
% right_circle_x = linspace(0, sqrt(b) - a, 500);
% 
% left_upper_circle_y = sqrt(b - (left_circle_x - a).^2);
% left_lower_circle_y = -left_upper_circle_y;
% 
% right_upper_circle_y = sqrt(b - (right_circle_x + a).^2);
% right_lower_circle_y = -right_upper_circle_y;
% 
% X0_x = [left_circle_x, right_circle_x, flip(right_circle_x), flip(left_circle_x)];
% X0_y = [left_upper_circle_y, right_upper_circle_y, flip(right_lower_circle_y), flip(left_lower_circle_y)];

% Область Х1

while true
%     figure('Position', [100 100 1000 1000]);
%     hold on
%     patch(X0_x, X0_y, 'g'); 
%     plot(x1(1), x1(2), 'r*', 'MarkerSize', 10);

    % Вычисляем начальное приближение psi_0
    
    alpha_n = linspace(alpha_start, alpha_end, alpha_n_point_cnt);
    psi_0(1, :) = cos(alpha_n);
    psi_0(2, :) = sin(alpha_n);

    % Задача в обратном времени: вычисление сопряжённой переменной psi для каждого psi_0

    opts = odeset('RelTol', eps, 'AbsTol', eps, 'MaxStep', ode45_max_dt_step);

    for psi_idx = 1 : size(psi_0, 2)
        sol = ode45(@conj_sys, [-T, -t0], psi_0(:, psi_idx), opts);
    
        psi(psi_idx, :, 1) = interp1(sol.x, sol.y(1, :), t_grid);
        psi(psi_idx, :, 2) = interp1(sol.x, sol.y(2, :), t_grid);
    end

    % Задача в обратном времени: вычисление оптимального управления из условия максимума

    if abs(det(D)) <= eps
        D = D + eye(2) * 1e-2;
    end

    for psi_0_idx = 1 : size(psi_0, 2)
        for t_idx = 1 : size(t_grid, 2)
           res = P_supp_vec(D' * [psi(psi_0_idx, t_idx, 1); psi(psi_0_idx, t_idx, 2)]);
           
           opt_control(psi_0_idx, t_idx, 1) = res(1); 
           opt_control(psi_0_idx, t_idx, 2) = res(2); 
        end
    end

    % Задача в обратном времени: вычисление траекторий

    opts = odeset('RelTol', eps, 'AbsTol', eps, 'MaxStep', ode45_max_dt_step, 'Events', @events_handler);
   
    for psi_0_idx = 1 : size(psi_0, 2)
        sol = ode45(@(t, y) traj_sys(t, y, t_grid, opt_control(psi_0_idx, :, :)), ...
            [-T, -t0], x1, opts);
    
        y(psi_0_idx, :, 1) = interp1(sol.x, sol.y(1, :), t_grid, 'spline');
        y(psi_0_idx, :, 2) = interp1(sol.x, sol.y(2, :), t_grid, 'spline');  
    
        % пересекла ли траектория множество Х0
        if numel(sol.xe) > 0 
            % отбираем в t_grid те времена, которые <= t_opt
            t_seg_idx = find(t_grid <= min(sol.xe));
            
            t_end = size(t_seg_idx, 2);
            if size(t_seg_idx, 2) < size(t_grid, 2)
                t_end = t_end + 1;
            end
                 
            t_grid(t_end) = sol.xe(1);
            y(psi_0_idx, t_end, 1) = sol.ye(1, 1);
            y(psi_0_idx, t_end, 2) = sol.ye(2, 1);
            
            if psi_0_idx_opt == 0 || min(sol.xe) <= t_opt
                t_opt = min(sol.xe);
                psi_0_idx_opt = psi_0_idx;
                t_opt_idx = t_end;
            end
            
%             plot(y(psi_0_idx, 1 : t_end, 1), y(psi_0_idx, 1 : t_end, 2), 'k');
%             plot(y(psi_0_idx, t_end, 1), y(psi_0_idx, t_end, 2), 'k', 'MarkerSize', 5);
%         else
%             t_end = size(t_grid, 2);
        end
    end
    
    if (psi_0_idx_opt > 0)
%         plot(y(psi_0_idx_opt, 1 : t_opt_idx, 1), y(psi_0_idx_opt, 1 : t_opt_idx, 2), 'r');
%         plot(y(psi_0_idx_opt, t_opt_idx, 1), y(psi_0_idx_opt, t_opt_idx, 2), 'r', 'MarkerSize', 10);
    else
        disp('Задача быстродействия неразрешима при заданных ограничениях.');
        is_solved = false;
        break;
    end
    
%     axis equal;
%     title('X_0, X_1 and optimal trajectory (red colour)');
%     xlabel('x_1');
%     ylabel('x_2');
%     legend({'X0', 'X1'});
%     hold off;
    % Погрешность выполнения условия трансверсальности на правом конце

    psi_opt = [psi(psi_0_idx_opt, t_opt_idx, 1), psi(psi_0_idx_opt, t_opt_idx, 2)];
    single_psi_opt = psi_opt / norm(psi_opt);

    trans_cond_abs_error = abs(X0_supp_func(-single_psi_opt) - ...
        scalar_prod(-single_psi_opt, y(psi_0_idx_opt, t_opt_idx, :)));
    
    if (trans_cond_abs_error > eps && alpha_end - alpha_start > alpha_seg_lim)
        alpha_opt = alpha_n(psi_0_idx_opt);
        alpha_seg_len = (alpha_end - alpha_start) / 2;
        alpha_start = alpha_opt - alpha_seg_len / 2;
        alpha_end = alpha_opt + alpha_seg_len / 2;
        
        t_opt = -T - 1;
        t_opt_idx = 0;
        psi_0_idx_opt = 0;
    else 
        break;
    end
    
end

if is_solved

% Результаты

disp('Решение линейной задачи быстродействия:');
disp(['Время быстродействия T = ', num2str(t_opt + T)]);
disp(['Точка прибытия в X_0: (', num2str(y(psi_0_idx_opt, t_opt_idx, 1)), ', ', num2str(y(psi_0_idx_opt, t_opt_idx, 2)), ')' ]);
disp(['Погрешность выполнения условия трансверсальности в точке прибытия: ', num2str(trans_cond_abs_error)]);

%% Строим графики

%Строим область Х0

left_circle_x = linspace(a - sqrt(b), 0, 500);
right_circle_x = linspace(0, sqrt(b) - a, 500);

left_upper_circle_y = sqrt(b - (left_circle_x - a).^2);
left_lower_circle_y = -left_upper_circle_y;

right_upper_circle_y = sqrt(b - (right_circle_x + a).^2);
right_lower_circle_y = -right_upper_circle_y;

X0_x = [left_circle_x, right_circle_x, flip(right_circle_x), flip(left_circle_x)];
X0_y = [left_upper_circle_y, right_upper_circle_y, flip(right_lower_circle_y), flip(left_lower_circle_y)];

figure('Position', [100 100 1000 1000]);

hold on

patch(X0_x, X0_y, 'g'); 

% Область Х1

plot(x1(1), x1(2), 'r*', 'MarkerSize', 10);

% Фазовый портрет оптимальной траектории

plot(y(psi_0_idx_opt, 1 : t_opt_idx, 1), y(psi_0_idx_opt, 1 : t_opt_idx, 2), 'r');
plot(y(psi_0_idx_opt, t_opt_idx, 1), y(psi_0_idx_opt, t_opt_idx, 2), 'r.', 'MarkerSize', 10);

axis equal;
title('X_0, X_1 and optimal trajectory');
xlabel('x_1');
ylabel('x_2');
legend({'X0', 'X1', 'optimal trajectory', 'arrival point'});
hold off;

% Область P

x = linspace(alpha - sqrt(r), alpha + sqrt(r), 1000);

y_upper = ((r - (x - alpha).^2) / 9).^0.25;
y_lower = -y_upper;

P_x = [x, flip(x)];
P_y = [y_upper, flip(y_lower)];

figure('Position', [100 100 1000 1000]);
hold on;
patch(P_x, P_y, 'g');

% Фазовый портрет оптимального управления

plot(opt_control(psi_0_idx_opt, (1 : t_opt_idx), 1), opt_control(psi_0_idx_opt, (1 : t_opt_idx), 2), 'r', 'LineWidth', 3);

title('P and optimal control');
xlabel('u_1');
ylabel('u_2');
legend({'P', 'optimal control'});
axis ([alpha - sqrt(r) - 1, alpha + sqrt(r) + 1, -(r / 9)^0.25 - 1, (r / 9)^0.25 + 1]);
hold off;

% Фазовый портрет сопряжённых переменных
figure('Position', [100 100 1000 1000]);
plot(psi(psi_0_idx_opt, (1 : t_opt_idx), 1), psi(psi_0_idx_opt, (1 : t_opt_idx), 2));
title('psi_2(psi_1)');
xlabel('psi_1');
ylabel('psi_2');

% Компоненты сопряжённых переменных
figure('Position', [100 100 1000 1000]);

subplot(2, 1, 1);
plot(1 : t_opt_idx, psi(psi_0_idx_opt, (1 : t_opt_idx), 1));
title('psi_1(t)');
xlabel('t');
ylabel('psi_1');

subplot(2, 1, 2);
plot(1 : t_opt_idx, psi(psi_0_idx_opt, (1 : t_opt_idx), 2));
title('psi_2(t)');
xlabel('t');
ylabel('psi_2');

% Компоненты оптимальной траектории

figure('Position', [100 100 1000 1000]);

subplot(2, 1, 1);
plot(1 : t_opt_idx, y(psi_0_idx_opt, (1 : t_opt_idx), 1));
title('x_1(t)');
xlabel('t');
ylabel('x_1');

subplot(2, 1, 2);
plot(1 : t_opt_idx, y(psi_0_idx_opt, (1 : t_opt_idx), 2));
title('x_2(t)');
xlabel('t');
ylabel('x_2');

% Компоненты оптимального управления

figure('Position', [100 100 1000 1000]);

subplot(2, 1, 1);
plot(1 : t_opt_idx, opt_control(psi_0_idx_opt, (1 : t_opt_idx), 1));
title('u_1(t)');
xlabel('t');
ylabel('u_1');

subplot(2, 1, 2);
plot(1 : t_opt_idx, opt_control(psi_0_idx_opt, (1 : t_opt_idx), 2));
title('u_2(t)');
xlabel('t');
ylabel('u_2');

%%
end
end

function dpsi_dt = conj_sys(t, psi) 
    load oc_lw1_temp.mat C;
    
    dpsi_dt = (-C') * psi;
end


function dy_dt = traj_sys(t, y, t_grid, opt_control) 
    load oc_lw1_temp.mat C D g;
    u = [interp1(t_grid, opt_control(1, :, 1), t, 'spline'); 
         interp1(t_grid, opt_control(1, :, 2), t, 'spline') 
        ];
    
    dy_dt = C * y + D * u + g;
end


function [value, isterminal, direction] = events_handler(t, y)
    load oc_lw1_temp.mat is_in_X0;
        
    value = not(is_in_X0(y(1), y(2)));
    isterminal = 1;
    direction = 0;
end


function [supp_func] = X0_supp_func (l)
    load oc_lw1_temp.mat a b;
    
    l = [l(1), l(2)];
  
    if b > a^2
        if l(1) > 0 && ...
             abs(l(2)) <= (b - a * (a - l(1))) / sqrt(b - a^2)
        
            supp_func = sqrt(b) * norm(l) - a * l(1);
        elseif l(1) < 0 && ...
            abs(l(2)) <= (b - a * (a + l(1))) / sqrt(b - a^2)
        
            supp_func = sqrt(b) * norm(l) + a * l(1);
        elseif l(2) > 0 && ...
                abs(l(1)) <= (l(2) * sqrt(b - a^2) - b) / a + a
            
            supp_func = l(2) * sqrt(b - a^2);
        elseif l(2) < 0 && ...
            abs(l(1)) <= (abs(l(2)) * sqrt(b - a^2) - b) / a + a
        
            supp_func = -l(2) * sqrt(b - a^2);
        end
    elseif b == a^2 
        supp_func = 0;
    else 
        error('X0_supp_func: X0 is an empty set because b < a^2');
    end
end


function [supp_vec] = P_supp_vec (l)
    load oc_lw1_temp.mat eps r alpha;

    if (abs(l(2)) > eps) 
        par = roots([(18 * l(1) / l(2))^2, 9, 0, -r]);
        par = real(par((abs(imag(par)) < eps) & (real(par) >= 0)));
        
        supp_vec = [ 18 * (l(1) / abs(l(2))) * par.^1.5 + alpha; ...
                    sqrt(par) * sign(l(2)) ];
    else
        supp_vec = [sqrt(r) * sign(l(1)) + alpha; 0];
    end
end     