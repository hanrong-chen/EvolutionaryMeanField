# I. Process patient sequences from fasta file

####################################################################################################
# Input:
# * sequences, numseq, lenseq (from parse_fasta_headers)
# * consseq

# Output:
# * binarysequences (lenseq x numseq binary array)
# * noncons_vs_time : mutations w.r.t. consensus sequence at each time
# * noninit_vs_time : mutations w.r.t. initial sequence at each time

# 6/21/15: given a consensus, convert strings in sequences array to binary strings  
# also output an array with the locations of the mutations in each sequence  
# converted to fn parse_aa2binary
# 6/14/16: output both noninit_vs_time and noncons_vs_time
####################################################################################################

function parse_aa2binary(sequences,lenseq,numseq,consseq)
    binarysequences=zeros(Int,lenseq,numseq)
    noncons_vs_time=Array[]
    noninit_vs_time=Array[]

    for i=1:numseq
    	noncons=Int[]
		mutations=Int[]

        for j=1:lenseq
            if sequences[i][j]!=consseq[j]
                binarysequences[j,i]=1
                push!(noncons,j)
            end

			if binarysequences[j,i]!=binarysequences[j,1]
	            push!(mutations,j)
			end
        end

        push!(noncons_vs_time,noncons)
		push!(noninit_vs_time,mutations)
    end
    
    return binarysequences, noncons_vs_time, noninit_vs_time
end


####################################################################################################
# Output:
# * statesequences : array of VECTORS containing state numbers at site i (given in statesVec)
# * noncons_vs_time : mutations w.r.t. consensus sequence at each time
# * noninit_vs_time : mutations w.r.t. initial sequence at each time

# 2/26/16: parse amino acids in sequences[:] to state numbers
# 3/6/16: output both noncons_vs_time & noninit_vs_time
# 3/7/16: changed both to arrays of Int vectors
####################################################################################################

function parse_aa2states(sequences,lenseq,numseq,statesVec)
	statesequences=Array[]
	noncons_vs_time=Array[]
    noninit_vs_time=Array[]

	for i=1:numseq
		noncons=Int[]
		mutations=Int[]

		push!(statesequences,zeros(Int,lenseq))

		for j=1:lenseq
			statenum=findin(statesVec[j],sequences[i][j])
			if length(statenum)==0 # amino acid/gap wasn't found, which means it got bundled into "X"
				statenum=findin(statesVec[j],"X")
			end
			statesequences[i][j]=statenum[1] # "[1]" because statenum is a 1x1 array

			if statesequences[i][j]!=numstatesVec[j]
                push!(noncons,j)
			end

			if statesequences[i][j]!=statesequences[1][j]
                push!(mutations,j)
			end
		end

		push!(noncons_vs_time,noncons)
		push!(noninit_vs_time,mutations)
	end

	return statesequences, noncons_vs_time, noninit_vs_time
end