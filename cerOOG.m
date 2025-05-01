function [OOG_sol, OOG_value] = cerOOG(A, E, Cp, Fp, Cr, Fr, del, eps, optset)
    % Input
    % A: system matrix A
    % E: matrix input E
    % Cp: matrix performance output 
    % Fp: matrix input from attack
    % Cr: matrix output monnitor
    % Fr: matrix input from attack to monitor
    % del: alarm threshold
    % eps: maximum attack energy
    % optset: settings for YALMIP
    % Output
    % OOG_sol: solution to OOG from solving optimization problem
    % OOG_value: value of OOG
    n = size(A, 1); % size of the certain part
    na = size(E, 2); % number of attack signals

    P = sdpvar(n, n, 'symmetric');
    g = sdpvar(1, 1);
    p = sdpvar(1, 1);
    sl = sdpvar(n + na, 1);
   
    F = [];
    F = [F, g >= 0, p >= 0, sl >= 0];
    % F = [F, P >= 0];
    lmi = [A'*P+P*A, P*E;
           E'*P, -p*eye(na)] + [Cp'; Fp'] * [Cp'; Fp']' - g * [Cr'; Fr'] * [Cr'; Fr']';
    F = [F, lmi <= diag(sl)];
    h = g * del + p * eps + ones(1, n + na) * sl * 1e6;
    OOG_sol = optimize(F, h, optset);
    OOG_value = value(g * del + p * eps);
end