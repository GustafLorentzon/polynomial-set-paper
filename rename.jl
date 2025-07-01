function rename(c)
    c_str=String(c)
    new_str=c_str;
    if (startswith(c_str,"Ba"))
        n=parse(Int,(c_str[3:end]))-1;
        new_str="a$n";
    elseif (startswith(c_str,"Bb"))
        n=parse(Int,(c_str[3:end]))-1;
        new_str="b$n";
    elseif (startswith(c_str,"B"))
        n=parse(Int,(c_str[2:end]))+1;
        new_str="Q$n";
    elseif c==:y
        new_str="c";
    end
    return Symbol(new_str)
end
function rename_to_new_notation(a::Tuple)
    return (rename(a[1]),a[2])
end

function rename_to_new_notation(g::Compgraph)

    # Rename coeff keys
    kk=keys(g.coeffs)
    for key in kk
        new_key=rename(key);
        coeffs=g.coeffs[key];
        delete!(g.coeffs,key);
        g.coeffs[new_key]=coeffs
    end

    # Rename parent keys, operations and parent vectors
    kk=keys(g.parents)
    for key in kk
        new_key=rename(key);
        new_parents=rename.(g.parents[key])        
        delete!(g.parents,key);
        g.parents[new_key]=new_parents

        oper=g.operations[key]
        delete!(g.operations,key);
        g.operations[new_key]=oper;
    end

    # Rename outputs
    g.outputs[:]=rename.(g.outputs);

  
end
