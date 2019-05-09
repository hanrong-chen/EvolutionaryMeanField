using Distributions
include("handle_epitopes.jl")
include("get_fitness_ept.jl")

################################################################################
# 6/5/16: find_hEff_t_v1
# 9/20/16: implement newept
# 9/21/16: turn off newept for an epitope if its sum_marginals>=1
# 11/1/16: return sum of b's for epitopes where sum_marginals<1
################################################################################

function find_hEff_t_v1(L,Jbeta,x_t;
						ept_on="off",ept_start=0,ept_end=0,ept_field=0,ept_init=0,
						newept=0)
	b_excl=0

	# compute RHS of mean-field equation
	# a. intrinsic landscape
	hEff_t_v1=diag(Jbeta).*(1 .- x_t) .+ Jbeta*x_t

	# 2b. epitope term, use ept_init
	if ept_on=="on"
		# 10/31/16: call a separate function to manage the effects of epitopes
		hEff_t_v1, b_excl=handle_epitopes_excl(x_t,ept_start,ept_end,ept_field,ept_init,
									   hEff_t_v1,excl=newept)
	end

	return hEff_t_v1, b_excl
end

################################################################################
# 6/30/16: run_DMF_inf_Ising from run_DMF_inf_Ising_v2
# 8/30/16: add ept_binary as input
# 8/31/16: change to t=0...tf
# 9/20/16: implement newept
# 10/3/16: rename to run_DMF_Ising_2; finitepop
# 10/5/16: include additional factor representing no. of offspring of consensus
# 10/13/16: "finitepop=1" means N is fixed, "finitepop=2" means N varies
# 10/15/16: abort when N(t)=0 or >N_0
# 11/1/16: multiply by exp(b_excl_t) (IX.78)
# 11/17/16: fix calculation of N(t) using mean fitness and sum of hEff*x
# 12/12/16: allow ept_field to be a num_ept x (tf+1) matrix
# 12/13/16: input N_stop & t_zone as well
# 180519: vectorize
# 181030: rename to run_EMF_Ising_3()
# 181101: convert to new T=AB version
# 181101: input vector of mutation rates instead

# Input:
# * L, muArray, tf, Jbeta, minit
# * ept_on, ept_start, ept_end, ept_field, ept_init
# * newept
# * finitepop, N_0
# * F_0

# Output:
# * hEff
# * m
# * NVec
################################################################################

function run_EMF_Ising_3(L,muArray,tf,Jbeta,minit;
						 ept_on="off",ept_start=0,ept_end=0,ept_field=0,ept_init=[],
						 finitepop=0,F_0=0,N_0=0,t_zone=0,N_stop_early=Inf,N_stop=Inf,newept=0)
	
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

	# 181101: handle muArray being an L x 2 array
	if size(muArray)==(L,2)
	elseif length(muArray)==L
		muArray=muArray.*ones(2)'
	elseif length(muArray)==1
		muArray=muArray*ones(L,2)
	else
		println("Error: muArray is $(size(muArray,1)) x $(size(muArray,2)) which is not $L x 2.")
	end

	# 6/3/16: construct hEff
	hEff=zeros(L,tf+1)

	##################################################
	# 5/27/16: initialize m
	##################################################
	m=zeros(L,tf+1)
	m[:,1]=minit
	
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

	##################################################
	# step 0
	# 6/3/16: initialize hEff
	# 6/5/16: use find_hEff_t_v1
	##################################################
	hEff[:,1],=find_hEff_t_v1(L,Jbeta,minit;
						ept_on="on",ept_start=ept_start,ept_end=ept_end,ept_field=ept_field[:,1],ept_init=ept_init,
						newept=newept)

	# 10/3/16: initialize NVec
	if finitepop==0
		NVec=zeros(tf+1)
		NVec[1]=N_0
	end
	##################################################
	if finitepop==1
		NVec=N_0*ones(Int,tf+1)
	elseif finitepop==2
		NVec=zeros(Int,tf+1)
		NVec[1]=N_0
	end

	# shoot forward in time, computing N[t], m[t], then hEff[t]
	for t=2:tf+1
		##################################################
		# step 1
		# 10/3/16: find NVec[t]
		##################################################
		# Ztilde_t=exp.(hEff[:,t-1]).*m[:,t-1] + (1 .- m[:,t-1]) # 181101: old T=BA version
		Ztilde_t=(1 .- muArray[:,1]).*(1 .- m[:,t-1]) .+ muArray[:,2].*m[:,t-1] +
				  exp.(hEff[:,t-1]) .* ( muArray[:,1].*(1 .- m[:,t-1]) + (1 .- muArray[:,2]).*m[:,t-1] ) # 181101: new T=AB version

		# 181030: put pre-computation of m[:,t] here because we need this (NOT m[:,t-1]) to find NVect[t]!
		# mdet_nomut_t=exp.(hEff[:,t-1]).*m[:,t-1] ./ Ztilde_t # 181101: old T=BA version
		# mdet_t=(1 .- muArray[:,2]).*mdet_nomut_t + muArray[:,1].*(1 .- mdet_nomut_t) # 181101: old T=BA version
		mdet_t=exp.(hEff[:,t-1]) .* ( muArray[:,1].*(1 .- m[:,t-1]) + (1 .- muArray[:,2]).*m[:,t-1] ) ./ Ztilde_t # 181101: new T=AB version

		if finitepop==2
			# 10/5/16: multiply by factor representing no. of offspring of consensus
			EHminusHtilde=F_0

			# 11/17/16: fix calculation of N(t) using mean fitness and sum of hEff*x
			mean_fitness,=get_fitness_ept(Jbeta,mdet_t,ept_start,ept_end,ept_init,ept_field[:,t]) # 181030: changed to mdet_t
			EHminusHtilde+=mean_fitness[1] # 180519: "[1]" since mean_fitness is a vector
			EHminusHtilde-=sum(hEff[:,t-1].*mdet_t) # 181030: changed to mdet_t

			Zapprox_t=prod(Ztilde_t)
			Zapprox_t*=exp(EHminusHtilde)

			# # 181101: include finitepop==0 case
			# if finitepop==0
			# 	NVec[t]=NVec[t-1]*Zapprox_t
			# end
			##################################################
			# 180519: draw from Poisson distribution instead
			NVec[t]=rand(Poisson(NVec[t-1]*Zapprox_t))

			# 10/15/16: abort calculation when N(t)=0 or 12/13/16: (N(t)>N_stop and t>t_zone)
			# 2/13/17: or when N(t)>1.5*Npeak_est at ANY time
			# 4/6/18: slight modification of stopping criterion. Added "t<=t_zone"
			# 180607: change to N_stop_early
			if NVec[t]==0 || (t<=t_zone && NVec[t]>N_stop_early) || (t>t_zone && NVec[t]>N_stop)
				break
			end
		end

		##################################################
		# step 2
		# find m[i,t]
		##################################################
		if finitepop==0
			m[:,t]=mdet_t
		else
			for i=1:L
				m[i,t]=rand(Binomial(NVec[t],mdet_t[i]))/NVec[t] # 180519
			end
		end

		##################################################
        # step 3
        # 6/5/16: find hEff[i,t] using find_hEff_t_v1
        ##################################################
		hEff[:,t],=find_hEff_t_v1(L,Jbeta,m[:,t];
						ept_on="on",ept_start=ept_start,ept_end=ept_end,ept_field=ept_field[:,t],ept_init=ept_init,
						newept=newept)
	end

	return hEff, m, NVec
end