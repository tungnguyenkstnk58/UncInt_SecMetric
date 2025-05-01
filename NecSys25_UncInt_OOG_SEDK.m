close all
clear
clc
%%
na = 2;
nyr = 2;
load SEDKGrid.mat;
A = A * 1e-2;
Ase = A(1:195, 1:195);
Adk = A(196:end, 196:end);
nse = length(Ase);
ndk = length(Adk);
A_dk_2_se = A(1:nse, nse+1:end);
A_se_2_dk = A(nse+1:end, 1:nse);
Tse = 2 * eye(nse);
Tdk = 2 * eye(ndk);
LGse = Tse + diag(Ase * ones(nse, 1)) + diag(A_dk_2_se * ones(ndk, 1)) - Ase;
LGdk = Tdk + diag(Adk * ones(ndk, 1)) + diag(A_se_2_dk * ones(nse, 1)) - Adk;
G = digraph(A');
G_concomp = conncomp(G, 'Type', 'strong', 'OutputForm', 'cell');

% System matrices 
B1 = [];
C2 = [];
set_dk_2_se = [];
for i = 1 : size(A_dk_2_se,2)
    if min(A_dk_2_se(:, i) == zeros(nse,1)) == 0
       set_dk_2_se = [set_dk_2_se, i]; 
       non_zero_ind = find(A_dk_2_se(:, i));
       for j = 1 : size(non_zero_ind, 1)
           tem = zeros(1, ndk);
           tem(i) = 1;
           C2 = [C2; tem];
           tem = zeros(nse, 1);
           tem(non_zero_ind(j)) = A_dk_2_se(non_zero_ind(j), i);
           B1 = [B1, tem];
       end
    end
end
B2 = [];
C1 = [];
set_se_2_dk = [];
for i = 1 : size(A_se_2_dk,2)
    if min(A_se_2_dk(:, i) == zeros(ndk,1)) == 0
        set_se_2_dk = [set_se_2_dk, i];
        non_zero_ind = find(A_se_2_dk(:, i));
        for j = 1 : size(non_zero_ind, 1)
            tem = zeros(ndk, 1);
            tem(non_zero_ind(j)) = 1;
            B2 = [B2, tem];
            tem = zeros(1, nse);
            tem(i) = A_se_2_dk(non_zero_ind(j), i);
            C1 = [C1; tem];
        end
    end
end
nco = size(B1, 2); % connection
B1 = [zeros(nse, nco); B1];
C2 = [C2, zeros(nco, ndk)];
C1 = [C1, zeros(nco, nse)];
B2 = [zeros(ndk, nco); B2];
A1 = [zeros(nse), eye(nse);
      - LGse, - Tse];
A2 = [zeros(ndk), eye(ndk);
      - LGdk, - Tdk];   
a = diag(ones(nse,1)); 
b = diag(ones(nse-1,1),1); b(nse, 1) = 1;
Cp = [a-b, zeros(nse, nse)];
Dp = zeros(nse, nco);
Fp = zeros(nse, na);
Dr = zeros(nyr, nco);
Fr = zeros(nyr, na);
D1 = zeros(nco, nco);
F1 = zeros(nco, na);

D2 = zeros(nco, nco);
Cr = zeros(nyr, nse*2);
Cr(1, nse+1) = 1; % choose monitor node
Cr(2, nse+1) = 1; % choose monitor node

%% Compute gamme_2 based on model
del = 1;
eps = 1e2;
optset = sdpsettings('verbose', 0, 'solver', 'mosek');
nG2 = ndk;
nG1 = nse;

[sol, g2] = cerOOG(A2, B2, C2, D2, zeros(1, nG2*2), zeros(1, nco), 0, 1, optset);

% Simulation values
exp = 1;
OOG_upp = zeros(exp, 1);
upp_time = zeros(exp, 1);
OOG_mid = zeros(exp, 1);
mid_time = zeros(exp, 1);
OOG_true = zeros(exp, 1);
true_time = zeros(exp, 1);

for expk = 1 : exp    
    Fx = [zeros(nse,na); Ase(:, randi([1 nse], na, 1))]; 
    %% Compute uncertain OOG
    [sol_upp, sol_v] = unOOG(A1, B1, Fx, Cp, Dp, Fp, Cr, Dr, Fr, C1, D1, F1, del, eps, g2, optset);
    if min(sol_upp.info(1:19) == 'Successfully solved')
        OOG_upp(expk) = sol_v;
        upp_time(expk) = sol_upp.solvertime;
    else
        sol_upp.info
    end
    
    %% Compute true OOG
    iD = (eye(nco + nco) - [zeros(nco), D1;D2, zeros(nco)]) \ eye(nco + nco);
    Ag = blkdiag(A1, A2) + [zeros(nG1*2, nco), B1;B2, zeros(nG2*2, nco)] * iD * blkdiag(C1, C2);
    Eg = [Fx; zeros(nG2*2, na)] + [zeros(nG1*2, nco), B1;B2, zeros(nG2*2, nco)] * iD * [F1; zeros(nco, na)];
    Crg = [Cr, zeros(nyr, nG2*2)] + [zeros(nyr, nco), Dr] * iD * blkdiag(C1, C2);
    Frg = Fr + [zeros(nyr, nco), Dr] * iD * [F1; zeros(nco, na)];
    Cpg = [Cp, zeros(nG1, nG2*2)] + [zeros(nG1, nco), Dp] * iD * blkdiag(C1, C2);
    Fpg = Fp + [zeros(nG1, nco), Dp] * iD * [F1; zeros(nco, na)];
    
    [sol_true, sol_v] = cerOOG(Ag, Eg, Cpg, Fpg, Crg, Frg, del, eps, optset);
    if min(sol_true.info(1:19) == 'Successfully solved')
        OOG_true(expk) = sol_v; %value(g * del + p * eps);
        true_time(expk) = sol_true.solvertime;
    else
        sol_true.info
    end

    %% Result collection
    [OOG_upp(expk), OOG_mid(expk), OOG_true(expk)]
    [upp_time(expk), mid_time(expk), true_time(expk)]

end









