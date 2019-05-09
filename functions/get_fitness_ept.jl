################################################################################
# 6/20/16: fitness function
# 8/31/16: include sequence[i]!=0 condition, simplify computation of fitness
################################################################################

function get_fitness(Jbeta,sequence)
    fitness=0
    for i=1:L
        if sequence[i]!=0
            temp=0
            temp+=Jbeta[i,i]
            for j=i+1:L
                temp+=Jbeta[i,j]*sequence[j]
            end
            fitness+=temp*sequence[i]
        end
    end
    return fitness
end

################################################################################
# 6/21/16: fitness function with epitopes
# 8/31/16: include sequence[i]!=0 condition, simplify computation of fitness,
#          output fitness_noept
# 11/17/16: account for ept_field being a vector
# 180515: "sequence" -> "marginals", marginals is L x tf+1 matrix
# 181101: remove tf as an argument
# 181207: change tf+1 -> num_seq since that's simply the number of sequences for which we compute the fitness
# 181207: vectorize sequences when computing fitnesses
################################################################################

function get_fitness_ept(Jbeta,marginals,proteptstart,proteptend,proteptbinary,ept_field)
    num_ept=length(proteptstart)
    num_seq=size(marginals,2) # number of sequences for which to find the fitness
    fitness=zeros(num_seq)

  #   for i=1:L
		# fitness+=Jbeta[i,i]*marginals[i,:]
		# for j=i+1:L
		# 	fitness+=Jbeta[i,j]*marginals[i,:].*marginals[j,:]
		# end
  #   end

    # 181207: vectorize sequences when computing fitnesses
    for s=1:num_seq
        fitness[s]=0.5*marginals[:,s]'*Jbeta*marginals[:,s] + dot(diag(Jbeta),marginals[:,s].*(1 .- 0.5*marginals[:,s]))
    end

    fitness_noept=fitness[:]
    
    # 12/12/16: handle ept_field being a num_ept x num_seq matrix
    if size(ept_field)==(num_ept,num_seq)
    elseif length(ept_field)==num_seq
        if num_seq==num_ept
            println("Warning: num_seq = num_ept. Will assume that ept_field is different for each sequence.")
        end
        ept_field=ept_field'.*ones(num_ept)
    elseif length(ept_field)==num_ept
        ept_field=ept_field.*ones(num_seq)'
    elseif length(ept_field)==1
        ept_field=ept_field*ones(num_ept,num_seq)
    else
        println("Error: ept_field is $(size(ept_field,1)) x $(size(ept_field,2)) which is not $num_ept x $num_seq.")
    end

    # 6/21/16: include epitopes
    # 8/8/16: modify to allow calculation of mean fitness
    # 180515: make term_in_ept a length-tf vector
    for k=1:length(proteptstart)
        term_in_ept=ones(num_seq)
        for i=proteptstart[k]:proteptend[k]
            if proteptbinary[k][i-proteptstart[k]+1]==0
                term_in_ept.*=1 .- marginals[i,:]
            else
                term_in_ept.*=marginals[i,:]
            end
        end
        fitness+=ept_field[k,:].*term_in_ept

        # # 181207: vectorize
        # term_in_ept=prod((1-marginals[proteptstart[k]:proteptend[k],:]).^((1-proteptbinary[k]).*ones(num_seq)').*
        #                  marginals[proteptstart[k]:proteptend[k],:].^(proteptbinary[k].*ones(num_seq)'),1)

        # # ### 181207: ASSUME EPITOPE TARGET SEQUENCES ARE ALL ZEROES ###
        # # term_in_ept=prod(1-marginals[proteptstart[k]:proteptend[k],:],1)
        # # println("Warning: all epitope target sequences assumed zeroes.")
        
        # fitness+=ept_field[k,:].*term_in_ept'
    end
    
    return fitness, fitness_noept
end