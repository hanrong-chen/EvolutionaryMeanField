####################################################################################################
# Input:
# * file in fasta format (i.e. headers starting with '>' followed by a sequence separated by newlines)

# Output:
# * headersall
# * sequencesall
# Both are string arrays containing all the headers and sequences in the fasta file.

# 6/20/15: download headers and sequences from fasta files from John  
# 6/21/15: converted to fn parse_fasta_file
# 181102: renamed
####################################################################################################

function read_fasta_file(filename)

    f=open(filename)

    # define variables for headers and sequences

    headersall=AbstractString[]
    sequencesall=AbstractString[]
    tempstring="" # builds up the current sequence, until the next '>'

    # parse first line

    line=chomp(readline(f)) # chomp removes a trailing newline
    if line[1]=='>'
        push!(headersall,line[2:end])
    else
        println("ERROR: first line doesn't start with '>'.")
    end

    # parse subsequent lines

    while !eof(f)
        line=chomp(readline(f)) # chomp removes a trailing newline
        if line[1]!='>'
            tempstring=tempstring*line # "*" means concatenate
        else
            # reached the next header
            push!(sequencesall,tempstring)
            push!(headersall,line[2:end])
            tempstring="" # reset string
        end
    end

    # push the last sequence in

    push!(sequencesall,tempstring)

    close(f)
    
    return headersall, sequencesall
end

####################################################################################################
# 181102: remove HXB2 header and sequence (split up from parse_fasta_headers())
####################################################################################################
function remove_hxb2_fasta(headersall,sequencesall)
    if findfirst("HXB2",headersall[1]) != nothing
        println("The first header contains \"HXB2\". Returning only subsequent headers and sequences.")
        return headersall[2:end], sequencesall[2:end]
    else
        println("The first header does not contain \"HXB2\". Returning all headers and sequences.")
        return headersall[:], sequencesall[:]
    end
end

####################################################################################################
# To edit this code:
# * Each header contains information separated by '.'. Currently, timepoint is in 3rd place,
#   and number of samples is in 4th place. If this is different, change code.
# * Also add new variables if I want to save more information.

# Input:
# * headers

# Output:
# * timepoints (Int array)
# * numsamples (Int array)

# 6/20/15: parse the content in the headers  
# 6/21/15: converted to fn parse_fasta_headers
# 181102: split up into two and renamed
####################################################################################################

function parse_patient_headers(headers)
    
    # define arrays to contain timepoints and number of samples
    
    numseq=length(headers)
    timepoints=zeros(Int64,numseq)
    numsamples=zeros(Int64,numseq)

    # populate arrays with timepoints (3rd) and number of samples (4th)
    
    for i=1:numseq
        thisheader=split(headers[i],'.')
        timepoints[i]=parse(Int,thisheader[3])
        numsamples[i]=parse(Int,thisheader[4])
    end
    
    return timepoints, numsamples
end

####################################################################################################
# 181102: parse fasta headers from multiple sequence alignment
####################################################################################################
function parse_MSA_headers(headers)
    
    numseq=length(headers)
    seqnames=String[]
    patient_codes=String[]
    patient_IDs=String[]

    # populate arrays with patient_codes (5th)
    
    for i=1:numseq
        thisheader=split(headers[i],'.')
        push!(seqnames,thisheader[4])
        push!(patient_codes,thisheader[5])
        push!(patient_IDs,thisheader[6])
    end
    
    return seqnames, patient_codes, patient_IDs
end