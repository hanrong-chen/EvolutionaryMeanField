# 7/8/17: estimate time to extinction assuming no mutation
function estimate_t_ext(tf,N_0,F_0,num_ept,ept_field)
    if length(ept_field)!=tf+1
        println("Error: ept_field doesn't have tf+1 entries.")
        return
    end

    t_ext_est=tf+1
    N_t=N_0

    for t=0:tf
        growth_t=F_0+num_ept*ept_field[t+1]
        N_t*=exp(growth_t)
        if N_t<1
            t_ext_est=t
            break
        end
    end

    return t_ext_est
end