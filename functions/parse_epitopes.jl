# II. Process epitopes of a patient from epitope file

####################################################################################################
## `fn` parse_epitope_file

# To edit this code:
# * define more variables if I want to save more information from the epitope file

# Input:
# * file containing epitope information

# Format of John's epitope files is, by column:
# 1. protein ("Gag", "Env", "Pol", etc.)
# 2. start site
# 3. end site
# 4. epitope
# 5. start time
# 6. escape time
# 7. escaped (0/1)
# 8. AGP escape (0/1)
# 9. entropy
# 10. %M 1st
# 11. ???

# Output:
# * eptinprotein (#1) (String array)
# * eptstart (#2) (Int array)
# * eptend (#3) (Int array)
# * ept (#4) (String array)
# * starttime (#5) (Int array)

# 6/21/15: download and parse epitope files from John  
# converted to fn parse_epitope_file
# 9/9/16: return start time too
####################################################################################################

function parse_epitope_file(filename)

    epitope_data=readdlm(filename)

    # define variables for protein, start site, end site

    eptinprotein=AbstractString[]
    eptstart=Int[]
    eptend=Int[]
    ept=AbstractString[]
    starttime=Int[]

    # save variables

    for i=1:size(epitope_data,1)
        push!(eptinprotein,epitope_data[i,1])
        push!(eptstart,epitope_data[i,2])
        push!(eptend,epitope_data[i,3])
        push!(ept,epitope_data[i,4])
        push!(starttime,epitope_data[i,5])
    end
    
    return eptinprotein, eptstart, eptend, ept, starttime
end