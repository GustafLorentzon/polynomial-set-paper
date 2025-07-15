using LinearAlgebra
@inline function exp_sas_6mult_1div(A)
    T=promote_type(eltype(A),Int64)
    A_copy=similar(A,T); copyto!(A_copy, A);
    return exp_sas_6mult_1div!(A_copy)
end

@inline function exp_sas_6mult_1div!(A)
    T=promote_type(eltype(A),Int64) # Make it work for many 'bigger' types (matrices and scalars)
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
    # Computation order: A2 A4 A6 Ua Va Ub3 Ub Uc U Vb3 Vb V Z X P
    # Computing A2 with operation: mult
    mul!(memslots2,memslots1,memslots1)
    # Computing A4 with operation: mult
    mul!(memslots3,memslots2,memslots2)
    # Computing A6 with operation: mult
    mul!(memslots4,memslots2,memslots3)
    # Computing Ua with operation: lincomb
    # Computing Ua = x*I+x*A2+x*A4+x*A6
    coeff1=32382376266240000
    coeff2=1187353796428800
    coeff3=10559470521600
    coeff4=33522128640
    memslots5 .= coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4
    mul!(memslots5, true, I*coeff1, true, true)
    # Computing Va with operation: lincomb
    # Computing Va = x*I+x*A2+x*A4+x*A6
    coeff1=64764752532480000
    coeff2=7771770303897600
    coeff3=129060195264000
    coeff4=670442572800
    memslots6 .= coeff2.*memslots2 .+ coeff3.*memslots3 .+ coeff4.*memslots4
    mul!(memslots6, true, I*coeff1, true, true)
    # Computing Ub3 with operation: lincomb
    # Computing Ub3 = x*A2+x*A4+x*A6
    coeff1=40840800
    coeff2=16380
    coeff3=1
    memslots7 .= coeff1.*memslots2 .+ coeff2.*memslots3 .+ coeff3.*memslots4
    # Computing Ub with operation: mult
    mul!(memslots8,memslots7,memslots4)
    # Deallocating Ub3 in slot 7
    # Computing Uc with operation: lincomb
    # Computing Uc = x*Ub+x*Ua
    coeff1=1
    coeff2=1
    # Smart lincomb recycle Ub
    memslots8 .= coeff1.*memslots8 .+ coeff2.*memslots5
    # Deallocating Ua in slot 5
    # Computing U with operation: mult
    mul!(memslots5,memslots1,memslots8)
    # Deallocating A in slot 1
    # Deallocating Uc in slot 8
    # Computing Vb3 with operation: lincomb
    # Computing Vb3 = x*A2+x*A4+x*A6
    coeff1=1323241920
    coeff2=960960
    coeff3=182
    # Smart lincomb recycle A2
    memslots2 .= coeff1.*memslots2 .+ coeff2.*memslots3 .+ coeff3.*memslots4
    # Deallocating A4 in slot 3
    # Computing Vb with operation: mult
    mul!(memslots1,memslots2,memslots4)
    # Deallocating Vb3 in slot 2
    # Deallocating A6 in slot 4
    # Computing V with operation: lincomb
    # Computing V = x*Vb+x*Va
    coeff1=1
    coeff2=1
    # Smart lincomb recycle Vb
    memslots1 .= coeff1.*memslots1 .+ coeff2.*memslots6
    # Deallocating Va in slot 6
    # Computing Z with operation: lincomb
    # Computing Z = x*V+x*U
    coeff1=1
    coeff2=-1
    memslots2 .= coeff1.*memslots1 .+ coeff2.*memslots5
    # Computing X with operation: lincomb
    # Computing X = x*V+x*U
    coeff1=1
    coeff2=1
    # Smart lincomb recycle V
    memslots1 .= coeff1.*memslots1 .+ coeff2.*memslots5
    # Deallocating U in slot 5
    # Computing P with operation: ldiv
    memslots3 .=memslots2\memslots1
    # Deallocating Z in slot 2
    # Deallocating X in slot 1
    return memslots3 # Returning P
end

