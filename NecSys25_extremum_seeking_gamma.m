close all
clear
clc
%%
nG2 = 10;
ny2 = 1;
ny1 = 1;
rng(1)
T2 = 0.8 * eye(nG2);
[LG2, AG2] = ConnectedGraph(nG2, 0);
A2 = [zeros(nG2), eye(nG2);
      - (LG2 + T2), -T2];
B2 = zeros(nG2*2, ny1);
for i = 1 : ny1
    idx = randi([1 nG2], 1, 1);
    B2(nG2+idx, i) = 10 + 0.1 * randn(1, 1);
end 
C2 = zeros(ny2, nG2*2);
for i = 1 : ny2
    idx = randi([1 nG2], 1, 1);
    C2(i, idx) = 1;
end
D2 = zeros(ny2, ny1);
sys2 = ss(A2, B2, C2, D2);
% figure(3)
% sigma(sys2)
[sv,wout] = sigma(sys2);
[~,idx] = max(sv(1,:));
wout(idx)

%% Compute gamma_2 based on model

optset = sdpsettings('verbose', 0);
[sol, g2] = cerOOG(A2, B2, C2, D2, zeros(1, nG2*2), zeros(1, ny1), 0, 1, optset);

%% Compute gamma_2 based on ESC
Ts = 5e-2; Tsim = 1800; Ns = Tsim/Ts;
u1 = randn(ny1, Ns);
fre_input = randn(ny1, Ns);
uG2 = randn(ny1, Ns-1);
u2 = zeros(ny2, Ns);
x2 = zeros(nG2*2, Ns);
g2_esc = zeros(1, Ns-1);
eta = zeros(1, Ns);
per = zeros(ny1, Ns);
xi = zeros(ny1, Ns);
th = zeros(ny1, Ns);

wp = 0.015; % perturbation frequency
wu = 5; % interaction frequency
wh = wp * 1.6; % HPF frequency
wl = wp * 0.4; % LPF frequency
ki = 0.12; % integral gain
ktrunc = 500/Ts; % k truncated signal norm
   
for k = 1 : Ns - 1
    u1(:, k) = 0.5; 
    for i = 1 : ny1
        per(i, k) = 0.01 * sin(wp * Ts * k + 0 * pi/4) + th(i, k);
    end    
    fre_input(:, k) = u1(:, k) + per(:, k);% * Ts * k;
    uG2(:, k) = 2 * sin(fre_input(:, k) * Ts * k);
    x2(:, k+1) = (eye(nG2*2) + Ts * A2) * x2(:, k) + Ts * B2 * uG2(:, k);
    u2(:, k) = C2 * x2(:, k) + D2 * uG2(:, k);
    if k <= inf %ktrunc
        g2_esc(:, k) = norm(u2(:, 1:k))^2 / norm(uG2(:, 1:k))^2;
    else
        g2_esc(:, k) = norm(u2(:, k-ktrunc:k))^2 / norm(uG2(:, k-ktrunc:k))^2;
    end
    eta(:, k+1) = eta(:, k) - Ts * wh * eta(:, k) + Ts * wh * g2_esc(:, k);
    xi(:, k+1) = xi(:, k) - Ts * wl * xi(:, k) + Ts * wl * (g2_esc(:, k) - eta(:, k)) * per(:, k);     
    th(:, k+1) = th(:, k) + Ts * ki * xi(:, k);    
end

disp([g2, max(max(sv))^2, mean(g2_esc(end/2:end))])

%% Figure
figure(1)
set(gcf,'position',[2000,800,600,200]);
plot(Ts:Ts:Tsim-Ts, ones(1,Ns-1)*g2,'--','linewidth',3);
hold on
plot(Ts:Ts:Tsim-Ts, g2_esc,':','LineWidth',4);
ylim([0, g2 * 1.2]);
grid;
legend('$\gamma_u$','$\tilde{\gamma}_u(t)$','Interpreter','latex','FontSize',16,'Location','best');
xlabel('Time (s)','FontSize',16);
set(gca,'FontSize',16,'fontWeight','bold');
yticks(0:1:3);
xlim([0 Tsim]);
% print(gcf,'figs/computegamma','-djpeg','-r300');
figure(2)
set(gcf,'position',[2000,300,600,200]);
plot(Ts:Ts:Tsim-Ts, ones(1,Ns-1)*wout(idx),'--','linewidth',3);
hold on
plot(Ts:Ts:Tsim-Ts,fre_input(1,1:end-1),':','linewidth',4); 
grid;
legend('$\omega_u^\star$','$\omega_u(t)$','Interpreter','latex','FontSize',16,'Location','best');
xlabel('Time (s)','FontSize',16);
xlim([0 Tsim]);
set(gca,'FontSize',16,'fontWeight','bold');
% print(gcf,'figs/computeomega','-djpeg','-r300');
figure(3)
set(gcf,'position',[1300,300,600,400]);
plot(Ts:Ts:Tsim-Ts,u2(1,1:end-1),'-','linewidth',3);
hold on
plot(Ts:Ts:Tsim-Ts,uG2(1,1:end),':','linewidth',3);





