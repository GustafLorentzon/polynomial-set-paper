function output = dummy(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 y
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Ba3 with operation: lincomb
    coeff1 = -8.738087791091056;
    coeff2 = 1.0;
    Ba3 = coeff1*A + coeff2*B2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 11.655854424500696;
    coeff2 = 1.0;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 28.625408287186385;
    coeff2 = 111.22660884010244;
    coeff3 = 1.0;
    Ba4 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = -2.456228171871179;
    coeff2 = 1.0;
    Bb4 = coeff1*A + coeff2*B2;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba5 with operation: lincomb
    coeff1 = 2.5533028987112494;
    coeff2 = -13.391542213243468;
    coeff3 = -0.8398466825433142;
    coeff4 = 1.0;
    Ba5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 6.135797749700576;
    coeff2 = -164.83784433201728;
    coeff3 = -2.341149959594229;
    coeff4 = 1.0;
    Bb5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 20.618513224080218;
    coeff2 = -75.11050539643026;
    coeff3 = -0.4478993051527408;
    coeff4 = -0.43805576778161404;
    coeff5 = 1.0;
    Ba6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = 71.64236185340361;
    coeff2 = -868.3450009304865;
    coeff3 = -24.01682973722021;
    coeff4 = 21.403873757976374;
    coeff5 = 1.0;
    Bb7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = -1.0399337031710059;
    coeff2 = -102.02482312600773;
    coeff3 = -1.670400209993799;
    coeff4 = 1.0;
    Bb6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = -31.121150540649744;
    coeff2 = -29.298773646949837;
    coeff3 = -44.4018687574287;
    coeff4 = 64.4454908890449;
    coeff5 = -2.561430924086462;
    coeff6 = 1.0;
    Ba7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing y with operation: lincomb
    coeff1 = 1.0;
    coeff2 = 13.0;
    coeff3 = 240.00030437375548;
    coeff4 = -34999.98336852076;
    coeff5 = 50010.09705940852;
    coeff6 = 39.858580653446985;
    coeff7 = 1271.5724716042835;
    coeff8 = 9.87735116411186;
    y = coeff1*I + coeff2*A + coeff3*B2 + coeff4*B3 + coeff5*B4 + coeff6*B5 + coeff7*B6 + coeff8*B7;
    output = y;
end

