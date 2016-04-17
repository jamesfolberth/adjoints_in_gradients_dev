function [] = fval_plotter()

addpath('export_fig/');

%fval();
%psnr();
ssim();

end


function [] = fval()

savefile = 'resolution_bior4.4_sym_fval.mat';

S = load(savefile);
fun_all_pinv = S.fun_all_pinv;
fun_all = S.fun_all;
fun_all_star = S.fun_all_star;
fstar = min(fun_all_star(:));

figure();
iters = 1:length(fun_all);
clf();
hold on;
semilogy(iters, fun_all-fstar, 'LineWidth', 3);
semilogy(iters, fun_all_pinv-fstar, '--', 'LineWidth', 3);
set(gca, 'YScale', 'log', 'FontSize', 16);
hold off;
legend('true adjoint', 'pinv approx');
xlabel('Iteration');
ylabel('Objective value - optimal value');

end


function [] = psnr()

savefile = 'resolution_bior4.4_sym_psnr.mat';

S = load(savefile, 'fun_all_pinv', 'fun_all');
fun_all_pinv = S.fun_all_pinv;
fun_all = S.fun_all;

iters = 1:length(fun_all_pinv);
figure()
clf();
hold on;
plot(iters, fun_all, 'LineWidth', 3);
plot(iters, fun_all_pinv, '--', 'LineWidth', 3);
hold off
h = legend('true adjoint', 'pinv approx', 'Location', 'SouthEast');
set(h, 'FontSize', 16);
set(gca, 'FontSize', 16);

xlabel('Iteration');
ylabel('PSNR');

end

function [] = ssim()

savefile = 'resolution_bior4.4_sym_ssim.mat';

S = load(savefile, 'fun_all_pinv', 'fun_all');
fun_all_pinv = S.fun_all_pinv;
fun_all = S.fun_all;

iters = 1:length(fun_all_pinv);
figure()
clf();
hold on;
plot(iters, fun_all, 'LineWidth', 3);
plot(iters, fun_all_pinv, '--', 'LineWidth', 3);
hold off
h = legend('true adjoint', 'pinv approx', 'Location', 'SouthEast');
set(h, 'FontSize', 16);
set(gca, 'FontSize', 16);

xlabel('Iteration');
ylabel('SSIM');

%set(gcf, 'Color', 'w');
%export_fig 'pinv_adjoint_SSIM_bior4.4_sym.pdf';

end
