using LinearAlgebra
@inline function onediv_deg20(A)
    T=promote_type(eltype(A),Float64)
    A_copy=similar(A,T); copyto!(A_copy, A);
    return onediv_deg20!(A_copy)
end

@inline function onediv_deg20!(A)
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
    coeff1=0.2
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
    coeff1=-0.50405
    coeff2=-5.4
    coeff3=1.0
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4
    # Computing Bb4 with operation: lincomb
    # Computing Bb4 = x*A+x*B2+x*B3
    coeff1=0.472025
    coeff2=1.0
    coeff3=1.0
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4
    # Computing B4 with operation: mult
    mul!(memslots7,memslots3,memslots6)
    # Deallocating Ba4 in slot 3
    # Deallocating Bb4 in slot 6
    # Computing Ba6 with operation: lincomb
    # Computing Ba6 = x*A+x*B2+x*B3+x*B4
    coeff1=-0.7121496786422238
    coeff2=-12.896589908051501
    coeff3=-1.108877049819538
    coeff4=1.0
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots7
    # Computing Bb5 with operation: lincomb
    # Computing Bb5 = x*A+x*B2+x*B3+x*B4
    coeff1=3.951383979196362
    coeff2=26.510412890655196
    coeff3=9.913366357049144
    coeff4=1.0
    memslots6 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots7
    # Computing B5 with operation: mult
    mul!(memslots8,memslots5,memslots6)
    # Deallocating Ba5 in slot 5
    # Deallocating Bb5 in slot 6
    # Computing Bb6 with operation: lincomb
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=4.505480164593282
    coeff2=47.556725115240035
    coeff3=119.4561861624968
    coeff4=11.091510692770393
    coeff5=1.0
    memslots5 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots4 .+ coeff4.*memslots7 .+ coeff5.*memslots8
    # Computing B6 with operation: mult
    mul!(memslots6,memslots3,memslots5)
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 5
    # Computing y with operation: lincomb
    # Computing y = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=1.0
    coeff2=1.0
    coeff3=-37.92990946324242
    coeff4=-1856.032587633356
    coeff5=-287.4755394774351
    coeff6=-21.157096699401635
    coeff7=1.0
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

