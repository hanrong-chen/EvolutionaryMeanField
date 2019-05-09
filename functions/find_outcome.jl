# 12/14/16: find t1 just before N(t)>N_stop && t>=t_zone, or else return tf
# 1/18/17: also return first time when N(t)=0 if extinction occurred, and also "escape", "die" or "" depending on outcome
# 1/22/17: change to first time when N(t)>N_stop if escape occurred
# 180607: change to N_stop_early

function find_outcome(t_zone,tf,NVec,N_stop; N_stop_early=1e10)
    outcome=""
    t1=tf+1
    for t=t_zone:tf # we start at t_zone because it may be that N(t)>N_stop at early times
        if NVec[t+1]==0
            t1=t # 1/18/17: first time when N(t)=0
            outcome="die"
            break
        elseif NVec[t+1]>N_stop
            t1=t # 1/22/17: first time when N(t)>N_stop
            outcome="escape"
            break
        end
    end

    # 180607: include check for exceeding N_stop_early at times before t_zone
    for t=0:t_zone-1
        if NVec[t+1]>N_stop_early
            t1=t
            outcome="escape"
            break
        end
    end

    return outcome, t1
end