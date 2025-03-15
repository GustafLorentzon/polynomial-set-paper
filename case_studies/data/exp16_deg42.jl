using LinearAlgebra
@inline function exp16_deg42(A)
    T=promote_type(eltype(A),ComplexF64)
    A_copy=similar(A,T); copyto!(A_copy, A);
    return exp16_deg42!(A_copy)
end

@inline function exp16_deg42!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=10
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    memslots8 = similar(A,T)
    memslots9 = similar(A,T)
    memslots10 = similar(A,T)
    # Assign precomputed nodes memslots
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 Ba8 Bb8 B8 y
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing Ba3 with operation: lincomb
    # Computing Ba3 = x*A+x*B2
    coeff1=1.0 + 0.0im
    coeff2=0.3552243434787143 + 0.31301865003529306im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing Bb3 with operation: lincomb
    # Computing Bb3 = x*A+x*B2
    coeff1=0.8360830164249938 + 1.51375982815941im
    coeff2=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots5,memslots3,memslots4)
    # Deallocating Ba3 in slot 3
    # Deallocating Bb3 in slot 4
    # Computing Ba4 with operation: lincomb
    # Computing Ba4 = x*A+x*B3
    coeff1=1.9878700319712954 + 2.3140766205659444im
    coeff2=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots5
    # Computing Bb4 with operation: lincomb
    # Computing Bb4 = x*A+x*B2
    coeff1=-3.3481818283799867 - 0.18781404816961916im
    coeff2=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B4 with operation: mult
    mul!(memslots6,memslots3,memslots4)
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 4
    # Computing Ba5 with operation: lincomb
    # Computing Ba5 = x*A+x*B2+x*B3+x*B4
    coeff1=-6.343965400071081 + 2.322828637545328im
    coeff2=1.0 + 0.0im
    coeff3=4.7235917363318425 + 0.3893233288877068im
    coeff4=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing Bb5 with operation: lincomb
    # Computing Bb5 = x*A+x*B2+x*B3
    coeff1=2.1949415868635036 + 2.0339267529547174im
    coeff2=1.6408178687320785 + 0.01243887307297313im
    coeff3=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5
    # Computing B5 with operation: mult
    mul!(memslots7,memslots3,memslots4)
    # Deallocating Ba5 in slot 3
    # Deallocating Bb5 in slot 4
    # Computing Ba6 with operation: lincomb
    # Computing Ba6 = x*A+x*B3+x*B4+x*B5
    coeff1=-14.240898251796072 + 3.670042231513456im
    coeff2=22.946602606317207 + 2.501748494503562im
    coeff3=1.9035275263080027 + 0.9773304760815128im
    coeff4=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots5 .+ coeff3.*memslots6 .+ coeff4.*memslots7
    # Computing Bb7 with operation: lincomb
    # Computing Bb7 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-1.0916975778241045 + 2.513481630026917im
    coeff2=3.4655419645706314 - 12.374255408172536im
    coeff3=25.859631711089065 - 1.9611449316761307im
    coeff4=1.4328957125708943 + 0.3517948071250101im
    coeff5=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing Bb6 with operation: lincomb
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4
    coeff1=1.842230367505875 + 3.2280755835474673im
    coeff2=7.355447021031304 + 5.342469876166199im
    coeff3=2.535462535574459 + 0.6067500151302706im
    coeff4=1.0 + 0.0im
    memslots8 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing B6 with operation: mult
    mul!(memslots9,memslots3,memslots8)
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 8
    # Computing Ba7 with operation: lincomb
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-56.440070975556615 + 0.8038703427701969im
    coeff2=-197.1421793323268 - 74.31908859859942im
    coeff3=-144.83452884676714 - 156.19616597452864im
    coeff4=-23.871745569154317 - 16.82457097482417im
    coeff5=-6.922978371736924 - 6.045993005201862im
    coeff6=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9
    # Computing B7 with operation: mult
    mul!(memslots8,memslots3,memslots4)
    # Deallocating Ba7 in slot 3
    # Deallocating Bb7 in slot 4
    # Computing Ba8 with operation: lincomb
    # Computing Ba8 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=-4.536202610923116 + 13.557328060523977im
    coeff2=-1979.3896551330106 + 1462.6053392454762im
    coeff3=564.0422556918274 - 434.50769564034044im
    coeff4=8.144305475508778 + 250.9970217238186im
    coeff5=28.701901758801203 + 6.764740732410251im
    coeff6=6.940779666409219 - 11.712228051306315im
    coeff7=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9 .+ coeff7.*memslots8
    # Computing Bb8 with operation: lincomb
    # Computing Bb8 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-22.304569550343828 + 14.982632784676683im
    coeff2=-151.42095847042825 - 49.07623727449899im
    coeff3=58.114992519629375 + 27.790767654075474im
    coeff4=-13.172306847687533 + 9.572804867625475im
    coeff5=1.0958333873628239 + 1.646026205769386im
    coeff6=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9
    # Computing B8 with operation: mult
    mul!(memslots10,memslots3,memslots4)
    # Deallocating Ba8 in slot 3
    # Deallocating Bb8 in slot 4
    # Computing y with operation: lincomb
    # Computing y = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=1.0 - 1.3937897627641574e-146im
    coeff2=16.0 + 2.9854139966181868e-145im
    coeff3=2.2538293597212923e6 - 1.6990034794703126e6im
    coeff4=463255.667959996 - 204462.24043615483im
    coeff5=129768.81807976031 - 283412.16330084566im
    coeff6=16558.498460634066 - 37303.882941481046im
    coeff7=-8107.967740845825 + 13211.941108730856im
    coeff8=-895.5940131936953 - 171.6495357115899im
    coeff9=92.31655506179625 + 50.985953848847636im
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots9 .+ coeff8.*memslots8 .+ coeff9.*memslots10
    mul!(memslots1, true, I*coeff1, true, true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 9
    # Deallocating B7 in slot 8
    # Deallocating B8 in slot 10
    return memslots1 # Returning y
end

