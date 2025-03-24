using LinearAlgebra
@inline function exp13_deg30(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); copyto!(A_copy, A);
    return exp13_deg30!(A_copy)
end

@inline function exp13_deg30!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
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
    # Computation order: B2 Ba3 Bb3 B3 Ba4 Bb4 B4 Ba5 Bb5 B5 Ba6 Bb7 Bb6 B6 Ba7 B7 y
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing Ba3 with operation: lincomb
    # Computing Ba3 = x*A+x*B2
    coeff1=-8.738087791091056
    coeff2=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing Bb3 with operation: lincomb
    # Computing Bb3 = x*A+x*B2
    coeff1=11.655854424500696
    coeff2=1.0
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots5,memslots3,memslots4)
    # Deallocating Ba3 in slot 3
    # Deallocating Bb3 in slot 4
    # Computing Ba4 with operation: lincomb
    # Computing Ba4 = x*A+x*B2+x*B3
    coeff1=28.625408287186385
    coeff2=111.22660884010244
    coeff3=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5
    # Computing Bb4 with operation: lincomb
    # Computing Bb4 = x*A+x*B2
    coeff1=-2.456228171871179
    coeff2=1.0
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B4 with operation: mult
    mul!(memslots6,memslots3,memslots4)
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 4
    # Computing Ba5 with operation: lincomb
    # Computing Ba5 = x*A+x*B2+x*B3+x*B4
    coeff1=2.5533028987112494
    coeff2=-13.391542213243468
    coeff3=-0.8398466825433142
    coeff4=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing Bb5 with operation: lincomb
    # Computing Bb5 = x*A+x*B2+x*B3+x*B4
    coeff1=6.135797749700576
    coeff2=-164.83784433201728
    coeff3=-2.341149959594229
    coeff4=1.0
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing B5 with operation: mult
    mul!(memslots7,memslots3,memslots4)
    # Deallocating Ba5 in slot 3
    # Deallocating Bb5 in slot 4
    # Computing Ba6 with operation: lincomb
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=20.618513224080218
    coeff2=-75.11050539643026
    coeff3=-0.4478993051527408
    coeff4=-0.43805576778161404
    coeff5=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing Bb7 with operation: lincomb
    # Computing Bb7 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=71.64236185340361
    coeff2=-868.3450009304865
    coeff3=-24.01682973722021
    coeff4=21.403873757976374
    coeff5=1.0
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing Bb6 with operation: lincomb
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4
    coeff1=-1.0399337031710059
    coeff2=-102.02482312600773
    coeff3=-1.670400209993799
    coeff4=1.0
    memslots8 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing B6 with operation: mult
    mul!(memslots9,memslots3,memslots8)
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 8
    # Computing Ba7 with operation: lincomb
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=-31.121150540649744
    coeff2=-29.298773646949837
    coeff3=-44.4018687574287
    coeff4=64.4454908890449
    coeff5=-2.561430924086462
    coeff6=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9
    # Computing B7 with operation: mult
    mul!(memslots8,memslots3,memslots4)
    # Deallocating Ba7 in slot 3
    # Deallocating Bb7 in slot 4
    # Computing y with operation: lincomb
    # Computing y = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=1.0
    coeff2=13.0
    coeff3=240.00030437375548
    coeff4=-34999.98336852076
    coeff5=50010.09705940852
    coeff6=39.858580653446985
    coeff7=1271.5724716042835
    coeff8=9.87735116411186
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots5 .+ coeff5.*memslots6 .+ coeff6.*memslots7 .+ coeff7.*memslots9 .+ coeff8.*memslots8
    mul!(memslots1, true, I*coeff1, true, true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 5
    # Deallocating B4 in slot 6
    # Deallocating B5 in slot 7
    # Deallocating B6 in slot 9
    # Deallocating B7 in slot 8
    return memslots1 # Returning y
end

