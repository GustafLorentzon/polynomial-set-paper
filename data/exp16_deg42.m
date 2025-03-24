function output = dummy(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 Ba8 Bb8 B8 y
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.5844400794380722 + 1i*-1.3959605079718456;
    coeff2 = 1.0 + 0.0im;
    Ba3 = coeff1*A + coeff2*B2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 0.8362985632994195 + 1i*1.5133442880778576;
    coeff2 = 1.0 + 0.0im;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 6.375652247856689 + 1i*0.8854685182993819;
    coeff2 = 1.0 + 0.0im;
    Ba4 = coeff1*A + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = -3.3481818283799867 + 1i*-0.18781404816961916;
    coeff2 = 1.0 + 0.0im;
    Bb4 = coeff1*A + coeff2*B2;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba5 with operation: lincomb
    coeff1 = -6.808601692935942 + 1i*12.544857079714737;
    coeff2 = 1.569740480375011 + 1i*-1.4144292127681592;
    coeff3 = 4.724523222553011 + 1i*0.38997858506695776;
    coeff4 = 1.0 + 0.0im;
    Ba5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 6.3189377982077914 + 1i*0.1597771243574064;
    coeff2 = 2.6186647244428443 + 1i*-2.2702663901368747;
    coeff3 = 1.0 + 0.0im;
    Bb5 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = 8.249004837244222 + 1i*65.08368685504342;
    coeff2 = 39.85357502273375 + 1i*-28.077354006661096;
    coeff3 = 4.3822084556779854 + 1i*-1.1052993512751348;
    coeff4 = 1.0 + 0.0im;
    Ba6 = coeff1*A + coeff2*B3 + coeff3*B4 + coeff4*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = 10.510805190812118 + 1i*6.242097832373527;
    coeff2 = -52.78387748400507 + 1i*-22.248348093966708;
    coeff3 = 38.23646224322802 + 1i*-39.21905642772004;
    coeff4 = 2.762936078783901 + 1i*-1.4393789763917824;
    coeff5 = 1.0 + 0.0im;
    Bb7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 7.426908143003605 + 1i*2.5429417754486807;
    coeff2 = 19.097383764544848 + 1i*-1.8267678978456707;
    coeff3 = 2.5365103349743428 + 1i*0.6076881686091713;
    coeff4 = 1.0 + 0.0im;
    Bb6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = 304.85341078550323 + 1i*435.7943373663025;
    coeff2 = 465.39549800179304 + 1i*1928.7742248563882;
    coeff3 = -772.5735928721261 + 1i*553.373346132297;
    coeff4 = -87.8771857159332 + 1i*96.12479896319797;
    coeff5 = -19.413016446774918 + 1i*0.08626434714601752;
    coeff6 = 1.0 + 0.0im;
    Ba7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = -88.09458998969606 + 1i*-594.4012503240502;
    coeff2 = 46185.73848288826 + 1i*-92443.74038237116;
    coeff3 = -13022.104899843343 + 1i*5569.324256735599;
    coeff4 = 1090.9032066106483 + 1i*-4877.961600583514;
    coeff5 = -99.09941154910298 + 1i*-259.5681289585798;
    coeff6 = -47.937421360569786 + 1i*-37.2904613151398;
    coeff7 = 1.0 + 0.0im;
    Ba8 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6 + coeff7*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = 234.8135383752222 + 1i*94.63946355367607;
    coeff2 = 419.7945638907328 + 1i*1438.293227983436;
    coeff3 = 155.55535494032847 + 1i*-241.5123648028718;
    coeff4 = 34.97166608032268 + 1i*63.68491910279126;
    coeff5 = 4.034990110990166 + 1i*1.0781891932596022;
    coeff6 = 1.0 + 0.0im;
    Bb8 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing y with operation: lincomb
    coeff1 = 1.0 + 1i*2.4829354080981247e-146;
    coeff2 = 16.0 + 1i*-5.149370397595671e-145;
    coeff3 = 2.249960953975171e6 + 1i*-1.700002399650597e6;
    coeff4 = 228695.45187978813 + 1i*72397.45882504294;
    coeff5 = 134809.0105693524 + 1i*-60016.94543945243;
    coeff6 = 8762.785342197853 + 1i*2630.2339792540383;
    coeff7 = -677.1635203234329 + 1i*-1499.3980874377685;
    coeff8 = 17.147247082874106 + 1i*13.291164285010847;
    coeff9 = 0.26629380073811015 + 1i*-9.297753385893747e-150;
    y = coeff1*I + coeff2*A + coeff3*B2 + coeff4*B3 + coeff5*B4 + coeff6*B5 + coeff7*B6 + coeff8*B7 + coeff9*B8;
    output = y;
end

