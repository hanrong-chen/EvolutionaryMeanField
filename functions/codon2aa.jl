# 181101: triplet codon to amino acid function
# 181102: got rid of all instances of "else" so that if codon contains 'X', '-', etc. then
#         codon2aa() would not return anything
function codon2aa(codon)
    if codon[1]=='T'
        if codon[2]=='T'
            if codon[3]=='T' || codon[3]=='C'
                aa='F'
            elseif codon[3]=='A' || codon[3]=='G'
                aa='L'
            end
        elseif codon[2]=='C'
            aa='S'
        elseif codon[2]=='A'
            if codon[3]=='T' || codon[3]=='C'
                aa='Y'
            elseif codon[3]=='A' || codon[3]=='G'
                aa='*'
            end
        elseif codon[2]=='G'
            if codon[3]=='T' || codon[3]=='C'
                aa='C'
            elseif codon[3]=='A'
                aa='*'
            elseif codon[3]=='G'
                aa='W'
            end
        end
    elseif codon[1]=='C'
        if codon[2]=='T'
            aa='L'
        elseif codon[2]=='C'
            aa='P'
        elseif codon[2]=='A'
            if codon[3]=='T' || codon[3]=='C'
                aa='H'
            elseif codon[3]=='A' || codon[3]=='G'
                aa='Q'
            end
        elseif codon[2]=='G'
            aa='R'
        end
    elseif codon[1]=='A'
        if codon[2]=='T'
            if codon[3]=='G'
                aa='M'
            elseif codon[3]=='T' || codon[3]=='C' || codon[3]=='A'
                aa='I'
            end
        elseif codon[2]=='C'
            aa='T'
        elseif codon[2]=='A'
            if codon[3]=='T' || codon[3]=='C'
                aa='N'
            elseif codon[3]=='A' || codon[3]=='G'
                aa='K'
            end
        elseif codon[2]=='G'
            if codon[3]=='T' || codon[3]=='C'
                aa='S'
            elseif codon[3]=='A' || codon[3]=='G'
                aa='R'
            end
        end
    elseif codon[1]=='G'
        if codon[2]=='T'
            aa='V'
        elseif codon[2]=='C'
            aa='A'
        elseif codon[2]=='A'
            if codon[3]=='T' || codon[3]=='C'
                aa='D'
            elseif codon[3]=='A' || codon[3]=='G'
                aa='E'
            end
        elseif codon[2]=='G'
            aa='G'
        end
    end
    
    return aa
end