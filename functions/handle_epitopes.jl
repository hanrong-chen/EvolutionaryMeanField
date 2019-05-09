################################################################################
# 10/31/16: call a separate function to manage the effects of epitopes
# 10/31/16: handle ept_field being a vector of epitopes (not time)
# 11/1/16: return sum of b's for epitopes where sum_marginals<1
################################################################################

################################################################################
# excl=0 : normal mean-field
# excl=1 : assumes mutations within an epitope are mutually exclusive, UNTIL
# sum_marginals>=1 for that epitope, after which normal mean-field is followed
# We are doing this for each epitope SEPARATELY.
################################################################################

function handle_epitopes_excl(x_t,ept_start,ept_end,ept_field,ept_init,hEff_t;excl=0)
	# 10/31/16: handle ept_field being a vector
	if length(ept_field)==1
		ept_field=ept_field*ones(length(ept_start))
	end

	# 11/1/16: sum of b's for epitopes where sum_marginals<1
	b_excl=0

	# 3/3/16: handle multiple epitopes (indexed by k)
	for k=1:length(ept_start)
		if excl!=0
			# 9/20/16: find sum of x_i' and (1-x_i'') for this epitope
			sum_marginals=0
			for i=ept_start[k]:ept_end[k]
				if ept_init[k][i-ept_start[k]+1]==0 # xinit[i]=consensus
					sum_marginals+=x_t[i]
				else # xinit[i]=mutant
					sum_marginals+=1-x_t[i]
				end
			end
		else
			sum_marginals=1
		end

		# 9/20/16: if sum_marginals<1, loop over sites in epitope,
		# and add or substract b depending on initial state
		if sum_marginals<1
			# 11/1/16: sum of b's for epitopes where sum_marginals<1
			b_excl+=ept_field[k]

			for i=ept_start[k]:ept_end[k]
				if ept_init[k][i-ept_start[k]+1]==0 # xinit[i]=consensus
					hEff_t[i]+=abs(ept_field[k])
				else # xinit[i]=mutant
					hEff_t[i]-=abs(ept_field[k])
				end
			end
		else # sum_marginals>=1; follow normal mean-field
			for i=ept_start[k]:ept_end[k]
				temp_ept_field=abs(ept_field[k])
				for j=ept_start[k]:ept_end[k]
					if j!=i
						if ept_init[k][j-ept_start[k]+1]==0 # xinit[j]=consensus
							temp_ept_field*=1-x_t[j]
						else # xinit[j]=mutant
							temp_ept_field*=x_t[j]
						end
					end
				end
				# 5/25/16: check ept_init at site i
				if ept_init[k][i-ept_start[k]+1]==0 # xinit[i]=consensus
					hEff_t[i]+=temp_ept_field
				else # xinit[i]=mutant
					hEff_t[i]-=temp_ept_field
				end
			end
		end
	end

	return hEff_t, b_excl
end