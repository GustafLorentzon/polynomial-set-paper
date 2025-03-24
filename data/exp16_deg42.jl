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
    coeff1=1.5844400794380722 - 1.3959605079718456im
    coeff2=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing Bb3 with operation: lincomb
    # Computing Bb3 = x*A+x*B2
    coeff1=0.8362985632994195 + 1.5133442880778576im
    coeff2=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2
    # Computing B3 with operation: mult
    mul!(memslots5,memslots3,memslots4)
    # Deallocating Ba3 in slot 3
    # Deallocating Bb3 in slot 4
    # Computing Ba4 with operation: lincomb
    # Computing Ba4 = x*A+x*B3
    coeff1=6.375652247856689 + 0.8854685182993819im
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
    coeff1=-6.808601692935942 + 12.544857079714737im
    coeff2=1.569740480375011 - 1.4144292127681592im
    coeff3=4.724523222553011 + 0.38997858506695776im
    coeff4=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing Bb5 with operation: lincomb
    # Computing Bb5 = x*A+x*B2+x*B3
    coeff1=6.3189377982077914 + 0.1597771243574064im
    coeff2=2.6186647244428443 - 2.2702663901368747im
    coeff3=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5
    # Computing B5 with operation: mult
    mul!(memslots7,memslots3,memslots4)
    # Deallocating Ba5 in slot 3
    # Deallocating Bb5 in slot 4
    # Computing Ba6 with operation: lincomb
    # Computing Ba6 = x*A+x*B3+x*B4+x*B5
    coeff1=8.249004837244222 + 65.08368685504342im
    coeff2=39.85357502273375 - 28.077354006661096im
    coeff3=4.3822084556779854 - 1.1052993512751348im
    coeff4=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots5 .+ coeff3.*memslots6 .+ coeff4.*memslots7
    # Computing Bb7 with operation: lincomb
    # Computing Bb7 = x*A+x*B2+x*B3+x*B4+x*B5
    coeff1=10.510805190812118 + 6.242097832373527im
    coeff2=-52.78387748400507 - 22.248348093966708im
    coeff3=38.23646224322802 - 39.21905642772004im
    coeff4=2.762936078783901 - 1.4393789763917824im
    coeff5=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7
    # Computing Bb6 with operation: lincomb
    # Computing Bb6 = x*A+x*B2+x*B3+x*B4
    coeff1=7.426908143003605 + 2.5429417754486807im
    coeff2=19.097383764544848 - 1.8267678978456707im
    coeff3=2.5365103349743428 + 0.6076881686091713im
    coeff4=1.0 + 0.0im
    memslots8 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6
    # Computing B6 with operation: mult
    mul!(memslots9,memslots3,memslots8)
    # Deallocating Ba6 in slot 3
    # Deallocating Bb6 in slot 8
    # Computing Ba7 with operation: lincomb
    # Computing Ba7 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=304.85341078550323 + 435.7943373663025im
    coeff2=465.39549800179304 + 1928.7742248563882im
    coeff3=-772.5735928721261 + 553.373346132297im
    coeff4=-87.8771857159332 + 96.12479896319797im
    coeff5=-19.413016446774918 + 0.08626434714601752im
    coeff6=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9
    # Computing B7 with operation: mult
    mul!(memslots8,memslots3,memslots4)
    # Deallocating Ba7 in slot 3
    # Deallocating Bb7 in slot 4
    # Computing Ba8 with operation: lincomb
    # Computing Ba8 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7
    coeff1=-88.09458998969606 - 594.4012503240502im
    coeff2=46185.73848288826 - 92443.74038237116im
    coeff3=-13022.104899843343 + 5569.324256735599im
    coeff4=1090.9032066106483 - 4877.961600583514im
    coeff5=-99.09941154910298 - 259.5681289585798im
    coeff6=-47.937421360569786 - 37.2904613151398im
    coeff7=1.0 + 0.0im
    memslots3 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9 .+ coeff7.*memslots8
    # Computing Bb8 with operation: lincomb
    # Computing Bb8 = x*A+x*B2+x*B3+x*B4+x*B5+x*B6
    coeff1=234.8135383752222 + 94.63946355367607im
    coeff2=419.7945638907328 + 1438.293227983436im
    coeff3=155.55535494032847 - 241.5123648028718im
    coeff4=34.97166608032268 + 63.68491910279126im
    coeff5=4.034990110990166 + 1.0781891932596022im
    coeff6=1.0 + 0.0im
    memslots4 .= coeff1.*memslots1 .+ coeff2.*memslots2 .+ coeff3.*memslots5 .+ coeff4.*memslots6 .+ coeff5.*memslots7 .+ coeff6.*memslots9
    # Computing B8 with operation: mult
    mul!(memslots10,memslots3,memslots4)
    # Deallocating Ba8 in slot 3
    # Deallocating Bb8 in slot 4
    # Computing y with operation: lincomb
    # Computing y = x*I+x*A+x*B2+x*B3+x*B4+x*B5+x*B6+x*B7+x*B8
    coeff1=1.0 + 2.4829354080981247e-146im
    coeff2=16.0 - 5.149370397595671e-145im
    coeff3=2.249960953975171e6 - 1.700002399650597e6im
    coeff4=228695.45187978813 + 72397.45882504294im
    coeff5=134809.0105693524 - 60016.94543945243im
    coeff6=8762.785342197853 + 2630.2339792540383im
    coeff7=-677.1635203234329 - 1499.3980874377685im
    coeff8=17.147247082874106 + 13.291164285010847im
    coeff9=0.26629380073811015 - 9.297753385893747e-150im
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

