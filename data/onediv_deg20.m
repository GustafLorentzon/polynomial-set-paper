function output = dummy(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B2 Bb3 B3 Ba4 Ba5 Bb4 B4 Ba6 Bb5 B5 Bb6 B6 y
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Bb3 with operation: lincomb
    coeff1 = 0.2;
    coeff2 = 1.0;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = B2 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 2.0;
    coeff2 = 1.0;
    Ba4 = coeff1*B2 + coeff2*B3;
    % Computing Ba5 with operation: lincomb
    coeff1 = -0.50405;
    coeff2 = -5.4;
    coeff3 = 1.0;
    Ba5 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 0.472025;
    coeff2 = 1.0;
    coeff3 = 1.0;
    Bb4 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba6 with operation: lincomb
    coeff1 = -0.7121496786422238;
    coeff2 = -12.896589908051501;
    coeff3 = -1.108877049819538;
    coeff4 = 1.0;
    Ba6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 3.951383979196362;
    coeff2 = 26.510412890655196;
    coeff3 = 9.913366357049144;
    coeff4 = 1.0;
    Bb5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 4.505480164593282;
    coeff2 = 47.556725115240035;
    coeff3 = 119.4561861624968;
    coeff4 = 11.091510692770393;
    coeff5 = 1.0;
    Bb6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing y with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 1.0;
    coeff3 = -37.92990946324242;
    coeff4 = -1856.032587633356;
    coeff5 = -287.4755394774351;
    coeff6 = -21.157096699401635;
    coeff7 = 1.0;
    y = coeff1*I + coeff2*A + coeff3*B2 + coeff4*B3 + coeff5*B4 + coeff6*B5 + coeff7*B6;
    output = y;
end

