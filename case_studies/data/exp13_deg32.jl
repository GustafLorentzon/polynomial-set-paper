using LinearAlgebra
@inline function exp13_deg32(A)
    T=promote_type(eltype(A),ComplexF64)
    A_copy=similar(A,T); copyto!(A_copy, A);
    return exp13_deg32!(A_copy)
end

@inline function exp13_deg32!(A)
    T=promote_type(eltype(A),ComplexF64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=9
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
    # Assign precomputed nodes memslots
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: B2 Ba3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 y
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing Ba3 with operation: lincomb
    # Computing Ba3 = x*A+x*B2
    coeff1=0.3076923076923077 + 0.0im
    coeff2=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots4,memslots3,memslots2)
    # Deallocating Ba3 in slot 3
    # Computing Ba4 with operation: lincomb
    # Computing Ba4 = x*A+x*B2+x*B3
    coeff1=0.6171352002562656 + 6.427188698809457im
    coeff2=-0.17427127240479315 + 0.8409922928038154im
    coeff3=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4
    # Computing Bb4 with operation: lincomb
    # Computing Bb4 = x*A+x*B3
    coeff1=0.28339261479945715 - 6.9779902898899im
    coeff2=1.0 + 0.0im
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots4
    # Computing B4 with operation: mult
    mul!(memslots6,memslots3,memslots5)
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 5
    # Computing Ba5 with operation: lincomb
    # Computing Ba5 = x*A+x*B2+x*B3+x*B4
    coeff1=-4.335339375086392 + 3.8638171271463415im
    coeff2=-49.61744744920391 + 4.395916760821174im
    coeff3=1.8058214927576657 - 0.4087004276644926im
    coeff4=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots6
    # Computing Bb5 with operation: lincomb
    # Computing Bb5 = x*A+x*B2+x*B3
    coeff1=1.139877371635806 + 0.8262023866206643im
    coeff2=1.8708743642284915 - 1.261488439205723im
    coeff3=1.0 + 0.0im
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4
    # Computing B5 with operation: mult
    mul!(memslots7,memslots3,memslots5)
    # Deallocating Ba5 in slot 3
    # Deallocating Bb5 in slot 5
    # Computing Ba6 with operation: lincomb
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-2.6063646537634724 + 3.9942905483291558im
    coeff2=-130.60164967302217 + 178.15336512579856im
    coeff3=10.224487991459235 - 17.646907247116918im
    coeff4=2.737197535249673 - 3.08488044817926im
    coeff5=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing Bb7 with operation: lincomb
    # Computing Bb7 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-16.862212773163503 - 1.509027994193272im
    coeff2=-69.20358895358133 + 121.13038070014684im
    coeff3=12.89653839358628 - 15.318528150224182im
    coeff4=1.6035914976917571 - 2.2394724626121305im
    coeff5=1.0 + 0.0im
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing Bb6 with operation: lincomb
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4
    coeff1=23.464747808401782 + 6.5626763978520914im
    coeff2=-15.727019402820218 - 5.172887208946824im
    coeff3=0.22145679809976196 + 0.9818986119105414im
    coeff4=1.0 + 0.0im
    memslots8 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots6
    # Computing B6 with operation: mult
    mul!(memslots9,memslots3,memslots8)
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 8
    # Computing Ba7 with operation: lincomb
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-2338.689418143228 - 287.3608995980036im
    coeff2=-8371.317685240483 + 6097.411850858722im
    coeff3=1158.3732361542977 - 1708.644389433381im
    coeff4=151.8556600500139 - 116.95309567735661im
    coeff5=19.630193173447562 + 30.29587810124931im
    coeff6=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9
    # Computing B7 with operation: mult
    mul!(memslots8,memslots3,memslots5)
    # Deallocating Ba7 in slot 3
    # Deallocating Bb7 in slot 5
    # Computing y with operation: lincomb
    # Computing y = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=1.0 + 0.0im
    coeff2=13.0 + 0.0im
    coeff3=-273458.2690503324 + 262170.0694913806im
    coeff4=29165.09704151681 - 160209.3409117471im
    coeff5=5543.173128579215 - 5950.1462642589895im
    coeff6=1411.0826641268723 + 1174.7948718770742im
    coeff7=171.69823035992385 - 7.436477849678896im
    coeff8=1.6827342204988953 + 0.0im
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots4 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots9 .+ coeff8.*memslots8
    mul!(memslots1, true, I*coeff1, true, true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 4
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 9
    # Deallocating B7 in slot 8
    return memslots1 # Returning y
end

