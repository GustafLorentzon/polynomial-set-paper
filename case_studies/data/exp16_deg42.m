function output = dummy(A)
    n = size(A,1);
    I = eye(n,n);
    % Computation order: B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 Ba8 Bb8 B8 y
    % Computing B2 with operation: mult
    B2 = A * A;
    % Computing Ba3 with operation: lincomb
    coeff1 = 1.0 + 0.0im;
    coeff2 = 0.3552243434787143 + 1i*0.31301865003529306;
    Ba3 = coeff1*A + coeff2*B2;
    % Computing Bb3 with operation: lincomb
    coeff1 = 0.8360830164249938 + 1i*1.51375982815941;
    coeff2 = 1.0 + 0.0im;
    Bb3 = coeff1*A + coeff2*B2;
    % Computing B3 with operation: mult
    B3 = Ba3 * Bb3;
    % Computing Ba4 with operation: lincomb
    coeff1 = 1.9878700319712954 + 1i*2.3140766205659444;
    coeff2 = 1.0 + 0.0im;
    Ba4 = coeff1*A + coeff2*B3;
    % Computing Bb4 with operation: lincomb
    coeff1 = -3.3481818283799867 + 1i*-0.18781404816961916;
    coeff2 = 1.0 + 0.0im;
    Bb4 = coeff1*A + coeff2*B2;
    % Computing B4 with operation: mult
    B4 = Ba4 * Bb4;
    % Computing Ba5 with operation: lincomb
    coeff1 = -6.343965400071081 + 1i*2.322828637545328;
    coeff2 = 1.0 + 0.0im;
    coeff3 = 4.7235917363318425 + 1i*0.3893233288877068;
    coeff4 = 1.0 + 0.0im;
    Ba5 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing Bb5 with operation: lincomb
    coeff1 = 2.1949415868635036 + 1i*2.0339267529547174;
    coeff2 = 1.6408178687320785 + 1i*0.01243887307297313;
    coeff3 = 1.0 + 0.0im;
    Bb5 = coeff1*A + coeff2*B2 + coeff3*B3;
    % Computing B5 with operation: mult
    B5 = Ba5 * Bb5;
    % Computing Ba6 with operation: lincomb
    coeff1 = -14.240898251796072 + 1i*3.670042231513456;
    coeff2 = 22.946602606317207 + 1i*2.501748494503562;
    coeff3 = 1.9035275263080027 + 1i*0.9773304760815128;
    coeff4 = 1.0 + 0.0im;
    Ba6 = coeff1*A + coeff2*B3 + coeff3*B4 + coeff4*B5;
    % Computing Bb7 with operation: lincomb
    coeff1 = -1.0916975778241045 + 1i*2.513481630026917;
    coeff2 = 3.4655419645706314 + 1i*-12.374255408172536;
    coeff3 = 25.859631711089065 + 1i*-1.9611449316761307;
    coeff4 = 1.4328957125708943 + 1i*0.3517948071250101;
    coeff5 = 1.0 + 0.0im;
    Bb7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5;
    % Computing Bb6 with operation: lincomb
    coeff1 = 1.842230367505875 + 1i*3.2280755835474673;
    coeff2 = 7.355447021031304 + 1i*5.342469876166199;
    coeff3 = 2.535462535574459 + 1i*0.6067500151302706;
    coeff4 = 1.0 + 0.0im;
    Bb6 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4;
    % Computing B6 with operation: mult
    B6 = Ba6 * Bb6;
    % Computing Ba7 with operation: lincomb
    coeff1 = -56.440070975556615 + 1i*0.8038703427701969;
    coeff2 = -197.1421793323268 + 1i*-74.31908859859942;
    coeff3 = -144.83452884676714 + 1i*-156.19616597452864;
    coeff4 = -23.871745569154317 + 1i*-16.82457097482417;
    coeff5 = -6.922978371736924 + 1i*-6.045993005201862;
    coeff6 = 1.0 + 0.0im;
    Ba7 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6;
    % Computing B7 with operation: mult
    B7 = Ba7 * Bb7;
    % Computing Ba8 with operation: lincomb
    coeff1 = -4.536202610923116 + 1i*13.557328060523977;
    coeff2 = -1979.3896551330106 + 1i*1462.6053392454762;
    coeff3 = 564.0422556918274 + 1i*-434.50769564034044;
    coeff4 = 8.144305475508778 + 1i*250.9970217238186;
    coeff5 = 28.701901758801203 + 1i*6.764740732410251;
    coeff6 = 6.940779666409219 + 1i*-11.712228051306315;
    coeff7 = 1.0 + 0.0im;
    Ba8 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6 + coeff7*B7;
    % Computing Bb8 with operation: lincomb
    coeff1 = -22.304569550343828 + 1i*14.982632784676683;
    coeff2 = -151.42095847042825 + 1i*-49.07623727449899;
    coeff3 = 58.114992519629375 + 1i*27.790767654075474;
    coeff4 = -13.172306847687533 + 1i*9.572804867625475;
    coeff5 = 1.0958333873628239 + 1i*1.646026205769386;
    coeff6 = 1.0 + 0.0im;
    Bb8 = coeff1*A + coeff2*B2 + coeff3*B3 + coeff4*B4 + coeff5*B5 + coeff6*B6;
    % Computing B8 with operation: mult
    B8 = Ba8 * Bb8;
    % Computing y with operation: lincomb
    coeff1 = 1.0 + 1i*-1.3937897627641574e-146;
    coeff2 = 16.0 + 1i*2.9854139966181868e-145;
    coeff3 = 2.2538293597212923e6 + 1i*-1.6990034794703126e6;
    coeff4 = 463255.667959996 + 1i*-204462.24043615483;
    coeff5 = 129768.81807976031 + 1i*-283412.16330084566;
    coeff6 = 16558.498460634066 + 1i*-37303.882941481046;
    coeff7 = -8107.967740845825 + 1i*13211.941108730856;
    coeff8 = -895.5940131936953 + 1i*-171.6495357115899;
    coeff9 = 92.31655506179625 + 1i*50.985953848847636;
    y = coeff1*I + coeff2*A + coeff3*B2 + coeff4*B3 + coeff5*B4 + coeff6*B5 + coeff7*B6 + coeff8*B7 + coeff9*B8;
    output = y;
end

