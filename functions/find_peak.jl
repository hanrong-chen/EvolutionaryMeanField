# 2/13/17: find tpeak_est, Npeak_est, t_zone, N_stop together
function find_peakest_stop(tf,N_0,F_0,num_ept,ept_field; N_thresh=1e3)
    if length(ept_field)!=tf+1
        println("Error: ept_field doesn't have tf+1 entries.")
    end

    tpeak_est=tf+1
    Npeak_est=0
    growth_exponent=0

    for t=0:tf
        growth_t=F_0+num_ept*ept_field[t+1]
        if growth_t>0
            growth_exponent+=growth_t
        else
            tpeak_est=t-1
            Npeak_est=N_0*exp(growth_exponent)
            break
        end
    end

    # 2/13/17: for the rest of the time, check for t_zone
    t_zone=0
    N_t=Npeak_est
    for t=tpeak_est+1:tf
        growth_t=F_0+num_ept*ept_field[t+1]
        N_t*=exp(growth_t)
        if N_t<N_thresh
            t_zone=t
            break
        end
    end

    N_stop=round(Int,N_thresh*10)

    return tpeak_est, Npeak_est, t_zone, N_stop
end

# 12/14/16: find time just before N(t) starts decreasing
# 1/18/17: also return N(tpeak)
function find_peak(t_zone,NVec)
    tpeak=0 # so that if I'm starting with large N_0 and constant b, this will return 0
    for t=0:t_zone
        if NVec[t+1]>NVec[t+2]
            tpeak=t
            break
        end
    end

    return tpeak, NVec[tpeak+1]
end

# # 12/15/16: estimate tpeak and Npeak
# # 1/18/17: rename to estimate_peak(), and move to find_peak.jl
# # ASSUMPTION: ept_field is a function of time, but is EQUAL for all epitopes
# function estimate_peak(tf,N_0,F_0,num_ept,ept_field)
#     if length(ept_field)!=tf+1
#         println("Error: ept_field doesn't have tf+1 entries.")
#     end
    
#     tpeak_est=tf+1
#     growth_exponent=0

#     for t=0:tf
#         growth_t=F_0+num_ept*ept_field[t+1]
#         if growth_t>0
#             growth_exponent+=growth_t
#         else
#             tpeak_est=t-1
#             break
#         end
#     end

#     return tpeak_est, N_0*exp(growth_exponent)
# end

# # 1/18/17: find t_zone, N_stop for pop2 simulations
# function find_stop(tpeak_est,Npeak_est)
#     t_zone=tpeak_est*2
#     N_stop=min(Npeak_est/10,5e4)
#     return t_zone, N_stop
# end