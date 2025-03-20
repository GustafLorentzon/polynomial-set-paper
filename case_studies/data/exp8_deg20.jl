using LinearAlgebra
@inline function exp8_deg20(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); copyto!(A_copy, A);
    return exp8_deg20!(A_copy)
end

@inline function exp8_deg20!(A)
    T=promote_type(eltype(A),Float64) # Make it work for many 'bigger' types (matrices and scalars)
    # max_memslots=8
    n=size(A,1)
    # The first slots are precomputed nodes [:A]
    memslots2 = similar(A,T)
    memslots3 = similar(A,T)
    memslots4 = similar(A,T)
    memslots5 = similar(A,T)
    memslots6 = similar(A,T)
    memslots7 = similar(A,T)
    memslots8 = similar(A,T)
    # Assign precomputed nodes memslots
    memslots1=A # overwrite A
    # Uniform scaling is exploited.
    # No matrix I explicitly allocated.
    # Computation order: B2 Bb3 B3 Ba4 Ba5 Bb4 B4 Ba6 Bb5 B5 Bb6 B6 y
    # Computing B2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing Bb3 with operation: lincomb
    # Computing Bb3 = x*A+x*B2
    coeff1=0.5
    coeff2=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots4,memslots2,memslots3)
    # Deallocating Bb3 in slot 3
    # Computing Ba4 with operation: lincomb
    # Computing Ba4 = x*B2+x*B3
    coeff1=2.0
    coeff2=1.0
    memslots3 .= coeff1.*memslots2 .+ coeff2.*memslots4
    # Computing Ba5 with operation: lincomb
    # Computing Ba5 = x*A+x*B2+x*B3
    coeff1=2.3374451754385963
    coeff2=-2.5625
    coeff3=1.0
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4
    # Computing Bb4 with operation: lincomb
    # Computing Bb4 = x*A+x*B2+x*B3
    coeff1=1.4484649122807018
    coeff2=1.0
    coeff3=1.0
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4
    # Computing B4 with operation: mult
    mul!(memslots7,memslots3,memslots6)
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 6
    # Computing Ba6 with operation: lincomb
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4
    coeff1=2.8309861554443847
    coeff2=8.7616118485412
    coeff3=5.123957592622475
    coeff4=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots7
    # Computing Bb5 with operation: lincomb
    # Computing Bb5 = x*A+x*B2+x*B3+x*B4
    coeff1=6.389966463262669
    coeff2=6.697361614351532
    coeff3=2.1451472591988834
    coeff4=1.0
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots7
    # Computing B5 with operation: mult
    mul!(memslots8,memslots5,memslots6)
    # Deallocating Ba5 in slot 5
    # Deallocating Bb5 in slot 6
    # Computing Bb6 with operation: lincomb
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=-2.458444697550387
    coeff2=-3.6724346694235876
    coeff3=16.044090747953085
    coeff4=7.557067023178642
    coeff5=1.0
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots7 .+ coeff5.*memslots8
    # Computing B6 with operation: mult
    mul!(memslots6,memslots3,memslots5)
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 5
    # Computing y with operation: lincomb
    # Computing y = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=1.0
    coeff2=8.0
    coeff3=-6.657689892163032
    coeff4=50.445902670306005
    coeff5=19.754729172913187
    coeff6=2.8090057997411706
    coeff7=0.47388735786811004
    # Smart lincomb recycle A
    memslots1 .= coeff2.*memslots1 .+ coeff3.*memslots2 .+ coeff4.*memslots4 .+ coeff5.*memslots7 .+ coeff6.*memslots8 .+ coeff7.*memslots6
    mul!(memslots1, true, I*coeff1, true, true)
    # Deallocating B2 in slot 2
    # Deallocating B3 in slot 4
    # Deallocating B4 in slot 7
    # Deallocating B5 in slot 8
    # Deallocating B6 in slot 6
    return memslots1 # Returning y
end

