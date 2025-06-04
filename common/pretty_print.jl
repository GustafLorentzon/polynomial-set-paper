
function pretty_print_row(a,x_symbol,mute)

    dof=0;
    for j=1:size(a,1)
        s=a[j]
        isthis_integer= try
            isinteger(s)
        catch
            false
        end

        if (isthis_integer && isreal(s) && real(s)>=0 && real(s)<=9)

            if !mute
                print(Int(s))
            end

        else
            if !mute
                print(x_symbol)
            end

            dof += 1
        end
        if !mute
            print(" ");
        end
    end
    return dof
end
function pretty_print(degopt;x_symbol="×",mute=false)
    (Ha,Hb,c)=get_degopt_coeffs(degopt)
    return pretty_print(Ha,Hb,c,x_symbol=x_symbol,mute=mute);
end

function pretty_print(Ha,Hb,c;x_symbol="×",mute=false)
    m=size(Ha,1);
    a_symbol=x_symbol;
    b_symbol=x_symbol;
    c_symbol=x_symbol;

    if (x_symbol == :abc)
        a_symbol="a";
        b_symbol="b";
        c_symbol="c";
    end

    dof=0;
    degvec=[1];
    for i=1:m;

        #@show Ha[i,1:i+1]
        dof1=pretty_print_row(Ha[i,1:i+1],a_symbol,mute);

        if !mute
            for j=1:2*(m-i)
                print(" ");
            end
            print("| ");
        end

        dof2=pretty_print_row(Hb[i,1:i+1],b_symbol,mute);

        dof += dof1 + dof2
        if !mute
        print(repeat(" ",2*m-2*i+1));
        end

        a_i=findlast(Ha[i,1:i+1] .!= 0);
        b_i=findlast(Hb[i,1:i+1] .!= 0);
        push!(degvec,degvec[a_i-1]+degvec[b_i-1]);

        if (!(mute))
            print(degvec[end])
            print("\n");
        end

    end
    if (!mute)
        println();
    end


    dofc=pretty_print_row(c,c_symbol,mute)
    dof += dofc;
    if (!mute)
        print("   nof $(x_symbol): $(dof)")
        println();
    end

    return (dof, degvec)
end


function latex_print_vals(degopt)
    (Ha,Hb,c)=get_degopt_coeffs(degopt);

    m=size(Ha,1);
    delim="&"
    for i=1:m
        for j=1:(m+1)
            if (!isinteger(Ha[i,j] ))
                println("\$a_{$i,$j}\$ $delim $(Ha[i,j]) \\\\");
            end
        end
    end
    for i=1:m
        for j=1:(m+1)
            if (!isinteger(Hb[i,j] ))
                println("\$b_{$i,$j}\$ $delim $(Hb[i,j]) \\\\");
            end
        end
    end
    for i=1:m+2
        if (!isinteger(c[i] ))
            println("\$c_{$i}\$ $delim $(c[i]) \\\\");
        end
    end
end



function polynomial_to_latex(p::Polynomial,v)
    coeffs_vals = coeffs(p)
    terms = []

    # Determine degrees
    degrees = degree(p):-1:0
    nonzero_terms = [(deg, coeffs_vals[deg + 1]) for deg in degrees if coeffs_vals[deg + 1] != 0]

    print("\\[")
    print("\\frac{\\partial p}{\\partial $v}=");
    
    if length(nonzero_terms) <= 4
        # Print everything if there are 4 or fewer terms
        terms_latex = join(["$(c)X^{$d}" for (d, c) in nonzero_terms], " + ")
        println("\\[ $terms_latex \\]")
        return
    end

    # Extract 4 highest-degree terms
    highest_terms = nonzero_terms[1:4]
    o_degree = highest_terms[end][1] - 1

    # Format terms
    for (d, c) in highest_terms
        sign = c < 0 ? "-" : "+"
        abs_c = abs(c)
        term = "$(abs_c)X^{$d}"
        push!(terms, (sign, term))
    end

    # Build LaTeX string
    latex_str = terms[1][2]
    for (sign, term) in terms[2:end]
        latex_str *= " $sign $term"
    end
    latex_str *= " + \\mathcal{O}(X^{$o_degree})"

    println("$latex_str \\]")
end
using Polynomials

function pretty_leading_terms(p::Polynomial)
    coeffs_vals = coeffs(p)
    degrees = degree(p):-1:0
    nonzero_terms = [(deg, coeffs_vals[deg + 1]) for deg in degrees if coeffs_vals[deg + 1] != 0]

    if isempty(nonzero_terms)
        return "0"
    end

    if length(nonzero_terms) <= 4
        return join(["$(c)X^$d" for (d, c) in nonzero_terms], " + ")
    end

    leading = nonzero_terms[1:4]
    o_degree = leading[end][1] - 1

    terms = ["$(c)X^$d" for (d, c) in leading]
    return join(terms, " + ") * " + O(X^$o_degree)"
end
