## @knitr ExtractIDs
ExtractIDs = function(probe_filter) {
    #probe_filter: a vector of S. cerevisiae genes
    #Get both S. pombe & S. cerevisiae ids from yeast2GENENAME library
    require(yeast2.db)
    genenames = as.list(yeast2GENENAME)
    probes = names(genenames) 
    
    #Get all transcript ids from yeast2annotation.csv
    annotations = read.csv(file='yeast2annotation.csv', header=TRUE, 
                           stringsAsFactors=FALSE)
    transcript_id = annotations[ ,3]
    probeset_id = annotations[ ,1]
    
    #Reorder the transcript_id to match probes
    transcript_id = transcript_id[match(probes, probeset_id)]
    
    #Retrieve the probeset and transcript ids for S. cerevisiae
    c_probe_id = probes[-match(probe_filter, probes)]
    c_transcript_id = transcript_id[-match(probe_filter, probes)]
    
    #We need the TranscriptID if the gene name is `NA'
    yeast_genenames = transcript_id
    for(i in seq(along=probeset_id)) {
        gname = genenames[i][[1]]
        if(!is.na(gname))
            yeast_genenames[i] = gname
    }
    
    #Set the gene name
    c_genename = yeast_genenames[-match(probe_filter, probes)]
    df = data.frame(probe=c_probe_id, transcript=c_transcript_id, 
                    genename=c_genename, stringsAsFactors=FALSE)
    return(df)
}