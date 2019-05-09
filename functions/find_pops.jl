function find_pops(L,muArray,tf,N_0,F_0,hEff_dmf1,m_dmf1,mean_fitness_vec)
    # 181101: handle muArray being an L x 2 array
    if size(muArray)==(L,2)
    elseif length(muArray)==L
        muArray=muArray.*ones(2)'
    elseif length(muArray)==1
        muArray=muArray*ones(L,2)
    else
        println("Error: muArray is $(size(muArray,1)) x $(size(muArray,2)) which is not $L x 2.")
    end
    
    # initialize pops
    pops=zeros(tf+1)
    pops[1]=N_0
    
    # construct Ztilde_matrix
    # Ztilde_matrix=exp.(hEff_dmf1).*m_dmf1 + (1-m_dmf1) # 181101: old T=BA version
    Ztilde_matrix=(1 .- muArray[:,1].*ones(tf+1)').*(1 .- m_dmf1) .+ (muArray[:,2].*ones(tf+1)').*m_dmf1 .+
                  exp.(hEff_dmf1) .* ( (muArray[:,1].*ones(tf+1)').*(1 .- m_dmf1) + (1 .- muArray[:,2].*ones(tf+1)').*m_dmf1 ) # 181101: new T=AB version
    
    log_Ztilde_matrix=log.(Ztilde_matrix)
    
    # compute last term on RHS of Gibbs inequality
    lastterm=zeros(tf+1)
    for t=2:tf+1
        lastterm[t]=sum(hEff_dmf1[:,t-1].*m_dmf1[:,t]) # see e.g. XV.143 comment
    end
        
    # find pops
    logZ=sum(log_Ztilde_matrix,dims=1)' .+ F_0 .+ mean_fitness_vec .- lastterm
    for t=2:tf+1
        pops[t]=pops[t-1]*exp(logZ[t])
    end
    
    return pops
end