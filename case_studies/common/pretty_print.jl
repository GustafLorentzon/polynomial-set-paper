
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
