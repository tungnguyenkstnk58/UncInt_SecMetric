function [OOG_sol, OOG_value] = unOOG(A1, B1, Fx, Cp, Dp, Fp, Cr, Dr, Fr, C1, D1, F1, del, eps, g2, optset)
    % Input
    % A1: system matrix A1
    % B1: matrix input B1 from uncertain part
    % Fx: matrix attack entry Fx
    % Cp: matrix output performance 
    % Dp: matrix input from uncertain part to performance
    % Fp: matrix input from attack
    % Cr: matrix output monnitor
    % Dr: matrix input from uncertain part to monitor
    % Fr: matrix input from attack to monitor
    % C1: matrix input for uncertain part from state of certain part
    % D1: matrix input for uncertain part from its output
    % F1: matrix input from attack to the input of uncertain part
    % del: alarm threshold
    % eps: maximum attack energy
    % g2: gain from uncertain part
    % optset: settings for YALMIP
    % Output
    % OOG_sol: solution to OOG from solving optimization problem
    % OOG_value: value of OOG
    n = size(A1, 1); % size of the certain part
    na = size(Fx, 2); % number of attack signals
    nco = size(B1, 2); % number of connections between uncertain and certain parts

    P = sdpvar(n, n, 'symmetric');
    g = sdpvar(1, 1);
    p = sdpvar(1, 1);
    t = sdpvar(1, 1);
    sl = sdpvar(n + na + nco, 1);
    
    F = [];
    F = [F, g >= 0, p >= 0, t >= 0, sl >= 0];
    % F = [F, P >= 0];
    lmi = [A1'*P+P*A1, P*Fx, P*B1;
           Fx'*P, -p*eye(na), zeros(na, nco);
           B1'*P, zeros(nco, na), -t*eye(nco)] + [Cp'; Fp'; Dp'] * [Cp'; Fp'; Dp']' ...
        - g * [Cr'; Fr'; Dr'] * [Cr'; Fr'; Dr']' ...
        + t * g2 * [C1'; F1'; D1'] * [C1'; F1'; D1']';
    F = [F, lmi <= diag(sl)];
    h = g * del + p * eps + ones(1, n + na + nco) * sl * 1e6;
    % h = h + 1e-2 * trace(P);
    OOG_sol = optimize(F, h, optset);
    OOG_value = value(g * del + p * eps);
end