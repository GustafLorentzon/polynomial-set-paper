function output = dummy(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B2 Ba3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 y
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Ba3 with operation: lincomb
    coeff1 = 0.3076923076923077 + 0.0im;
    coeff2 = 1.0 + 0.0im;
    Ba3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * B2;
    % Computing Ba4 with operation: lincomb
    coeff1 = 0.6171352002562656 + 1i*6.427188698809457;
    coeff2 = -0.17427127240479315 + 1i*0.8409922928038154;
    coeff3 = 1.0 + 0.0im;
    Ba4 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = 0.28339261479945715 + 1i*-6.9779902898899;
    coeff2 = 1.0 + 0.0im;
    Bb4 = coeff1*A + coeff2*B3;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba5 with operation: lincomb
    coeff1 = -4.335339375086392 + 1i*3.8638171271463415;
    coeff2 = -49.61744744920391 + 1i*4.395916760821174;
    coeff3 = 1.8058214927576657 + 1i*-0.4087004276644926;
    coeff4 = 1.0 + 0.0im;
    Ba5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 1.139877371635806 + 1i*0.8262023866206643;
    coeff2 = 1.8708743642284915 + 1i*-1.261488439205723;
    coeff3 = 1.0 + 0.0im;
    Bb5 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = -2.6063646537634724 + 1i*3.9942905483291558;
    coeff2 = -130.60164967302217 + 1i*178.15336512579856;
    coeff3 = 10.224487991459235 + 1i*-17.646907247116918;
    coeff4 = 2.737197535249673 + 1i*-3.08488044817926;
    coeff5 = 1.0 + 0.0im;
    Ba6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = -16.862212773163503 + 1i*-1.509027994193272;
    coeff2 = -69.20358895358133 + 1i*121.13038070014684;
    coeff3 = 12.89653839358628 + 1i*-15.318528150224182;
    coeff4 = 1.6035914976917571 + 1i*-2.2394724626121305;
    coeff5 = 1.0 + 0.0im;
    Bb7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 23.464747808401782 + 1i*6.5626763978520914;
    coeff2 = -15.727019402820218 + 1i*-5.172887208946824;
    coeff3 = 0.22145679809976196 + 1i*0.9818986119105414;
    coeff4 = 1.0 + 0.0im;
    Bb6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = -2338.689418143228 + 1i*-287.3608995980036;
    coeff2 = -8371.317685240483 + 1i*6097.411850858722;
    coeff3 = 1158.3732361542977 + 1i*-1708.644389433381;
    coeff4 = 151.8556600500139 + 1i*-116.95309567735661;
    coeff5 = 19.630193173447562 + 1i*30.29587810124931;
    coeff6 = 1.0 + 0.0im;
    Ba7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing y with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = 13.0 + 0.0im;
    coeff3 = -273458.2690503324 + 1i*262170.0694913806;
    coeff4 = 29165.09704151681 + 1i*-160209.3409117471;
    coeff5 = 5543.173128579215 + 1i*-5950.1462642589895;
    coeff6 = 1411.0826641268723 + 1i*1174.7948718770742;
    coeff7 = 171.69823035992385 + 1i*-7.436477849678896;
    coeff8 = 1.6827342204988953 + 0.0im;
    y = coeff1*I + coeff2*A + coeff3*B2 + coeff4*B3 + coeff5*B4 + coeff6*B5 + coeff7*B6 + coeff8*B7;
    output = y;
end

