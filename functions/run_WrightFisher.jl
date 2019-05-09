using Distributions
using StatsBase
include("get_fitness_ept.jl")

################################################################################
# 181207: constrain there to be at most one mutation per sequence

# Input:
# * seqs_init : Array of L x NVec[t] matrices
################################################################################
function run_WrightFisher(Jbeta,L::Int,mu,N::Int,rho,seqs_init,tf::Int;
                  ept_on="off",ept_start=0,ept_end=0,ept_field=0,ept_init=[])
    num_ept=length(ept_start)

    # 12/12/16: handle ept_field being a num_ept x (tf+1) matrix
    if size(ept_field)==(num_ept,tf+1)
    elseif length(ept_field)==tf+1
        if tf+1==num_ept
            println("Warning: tf+1 = num_ept. Will assume that ept_field is time-varying.")
        end
        ept_field=ept_field'.*ones(num_ept)
    elseif length(ept_field)==num_ept
        ept_field=ept_field.*ones(tf+1)'
    elseif length(ept_field)==1
        ept_field=ept_field*ones(num_ept,tf+1)
    else
        println("Error: ept_field is $(size(ept_field,1)) x $(size(ept_field,2)) which is not $num_ept x $(tf+1).")
    end

    # 180607: define minit for purposes of defining epitope sequence
    minit=sum(seqs_init,dims=2)/N

    if ept_on=="on"
        if ept_init==[]
            # 5/25/16: epitopes are the initial states. Construct ept_init
            ept_init=Array[] # contains vectors, each for one epitope

            # 3/3/16: handle multiple epitopes (indexed by k)
            for k=1:num_ept
                # ept_k=zeros(ept_end[k]-ept_start[k]+1)
                
                # 5/23/16: input 0 if minit[i]=consensus, 1 if minit[i]=mutant
                ept_k=Int.(minit[ept_start[k]:ept_end[k]].>=0.5) # mutant is at least as frequent initially
                
                push!(ept_init,ept_k)
            end
        else
            # # 8/30/16: use ept_binary if it is given
            # ept_init=ept_binary
        end
    end

    # 181206: initialize seqsMat
    seqsMat=zeros(L,N,tf+1)
    seqsMat[:,:,1]=seqs_init

    for t=2:tf+1
        # 2/4/17: calculate fitnesses of all sequences at time t-1
        fitnesses_parent,=get_fitness_ept(Jbeta,seqsMat[:,:,t-1],ept_start,ept_end,ept_init,ept_field[:,t-1])

        # 180607: correct computation of total number of offspring
        all_offspring=exp.(fitnesses_parent)
        sum_all_offspring=sum(all_offspring)

        # 180607: use StatsBase
        probabilities_all=all_offspring/sum(all_offspring)
        weights_all=Weights(probabilities_all)

        # 2/4/17: populate time t by drawing from seqsMat[:,:,t-1] weighted by weights_all
        seqs_t=zeros(L,N)
        for n=1:N
            seqs_t[:,n]=seqsMat[:,sample(collect(1:N),weights_all),t-1]

            # 181207: constrain there to be at most one mutation per sequence
            if rand()<mu*L
                s=rand(1:L)
                seqs_t[s,n]=1-seqs_t[s,n] # flip s
            end
        end

        # # 180607: implement mutations (OLD)
        # seqs_t=mod.(seqs_t + Int.(rand(size(seqs_t)).<mu),2)

        # 181207: implement recombinations
        num_crossings=rand(Poisson(rho*L*N))
        for c=1:num_crossings
            n1=rand(1:N)
            n2=rand(1:N)
            s1=rand(1:L-1)
            seqtemp1=seqs_t[:,n1]
            seqs_t[:,n1]=cat(seqs_t[1:s1,n1],seqs_t[s1+1:L,n2],dims=1)
            seqs_t[:,n2]=cat(seqs_t[1:s1,n2],seqtemp1[s1+1:L],dims=1)
        end

        # 181206: put into seqsMat
        seqsMat[:,:,t]=seqs_t
    end

    return seqsMat
end