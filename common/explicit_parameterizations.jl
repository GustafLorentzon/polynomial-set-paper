using LinearAlgebra

# Following function returns Triplet to compute a given 12 degree polynomial
function P12_coeff_solver(mc)
    c₆  = mc[12+1];
    a₃₃ = mc[11+1]/(2.0*c₆);
    a₃₂ = mc[10+1]/(2.0*c₆) - a₃₃^2/2.0 ;
    a₄₄ = mc[9+1]/(2.0*c₆) - 1/2.0 - a₃₂*a₃₃
    β₄₃ = mc[8+1]/c₆ - (a₃₃ + 2.0*a₃₃*a₄₄ + a₃₂^2)
    β₄₂ = mc[7+1]/c₆ - (a₃₂ + a₃₃*β₄₃ + 2.0*a₃₂*a₄₄);
    c₅  = mc[6+1]    - c₆*(a₄₄ + a₄₄^2 + a₃₃*β₄₂ + a₃₂*β₄₃);
    a₄₃ = mc[5+1]/c₆ - (a₃₃*c₅/c₆ + a₄₄*β₄₃ + a₃₂*β₄₂);
    a₄₂ = mc[4+1]/c₆ - (a₃₂*c₅/c₆ + a₄₄*β₄₂ + a₄₃*β₄₃ - a₄₃^2);
    c₄  = mc[3+1]    - c₆*(a₄₃*β₄₂ + a₄₂*β₄₃ - 2.0*a₄₂*a₄₃);
    c₃  = mc[2+1]    - c₆*(a₄₂*β₄₂ - a₄₂^2);
    c₂  = mc[1+1];
    c₁  = mc[0+1];

    b₄₂ = β₄₂ - a₄₂
    b₄₃ = β₄₃ - a₄₃
    

    HA = [0 1   0   0   0;
          0 1   0   0   0;
          0 a₃₂ a₃₃ 1   0;
          0 a₄₂ a₄₃ a₄₄ 1]

    HB = [0 1   0   0     0;
          0 0   1   0     0;
          0 0   0   1     0;
          0 b₄₂ b₄₃ a₄₄+1 1]

    C  = [c₁, c₂, c₃, c₄, c₅, c₆];
    return HA, HB, C
end





function P8_coeff_solver(mc)
      c₅  = mc[8+1]
      a₂₂ = mc[7+1] / (2*c₅)
      a₃₃ = mc[6+1] / (2*c₅) - (1 + a₂₂^2)/2
      β₃₂ = mc[5+1]/c₅ - (a₂₂ + 2*a₂₂*a₃₃)
      c₄  = mc[4+1] - c₅*(a₃₃ + a₃₃^2 + a₂₂*β₃₂)
      a₃₂ = mc[3+1]/c₅ - (a₂₂*c₄/c₅ + a₃₃*β₃₂)
      c₃  = mc[2+1] - c₅*(a₃₂*β₃₂ - a₃₂^2)
      c₂  = mc[1+1]
      c₁  = mc[0+1]

      b₃₂ = β₃₂ - a₃₂

      HA= [ 0 1   0   0;
            0 a₂₂ 1   0;
            0 a₃₂ a₃₃ 1];

      HB= [ 0 1   0     0;
            0 0   1     0;
            0 b₃₂ a₃₃+1 1]

      c = [c₁,c₂,c₃,c₄,c₅];
    return HA, HB, c
end