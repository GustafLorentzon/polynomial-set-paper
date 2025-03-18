function output = dummy(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B2 Bb3 B3 Ba4 Ba5 Bb4 B4 Ba6 Bb5 B5 Bb6 B6 y
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 0.5;
    coeff2 = 1.0;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = B2 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 2.0;
    coeff2 = 1.0;
    Ba4 = coeff1*B2 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = 2.3374451754385963;
    coeff2 = -2.5625;
    coeff3 = 1.0;
    Ba5 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 1.4484649122807018;
    coeff2 = 1.0;
    coeff3 = 1.0;
    Bb4 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba6 with operation: lincomb
    coeff1 = 2.8309861554443847;
    coeff2 = 8.7616118485412;
    coeff3 = 5.123957592622475;
    coeff4 = 1.0;
    Ba6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 6.389966463262669;
    coeff2 = 6.697361614351532;
    coeff3 = 2.1451472591988834;
    coeff4 = 1.0;
    Bb5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Bb6 with operation: lincomb
    coeff1 = -2.458444697550387;
    coeff2 = -3.6724346694235876;
    coeff3 = 16.044090747953085;
    coeff4 = 7.557067023178642;
    coeff5 = 1.0;
    Bb6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing y with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 8.0;
    coeff3 = -6.657689892163032;
    coeff4 = 50.445902670306005;
    coeff5 = 19.754729172913187;
    coeff6 = 2.8090057997411706;
    coeff7 = 0.47388735786811004;
    y = coeff1*I + coeff2*A + coeff3*B2 + coeff4*B3 + coeff5*B4 + coeff6*B5 + coeff7*B6;
    output = y;
end

