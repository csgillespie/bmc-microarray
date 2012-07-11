## @knitr RemoveProbes
RemoveProbes=function(listOutProbeSets, cdfpackagename, probepackagename){    
    #listOutProbeSets: Probes sets that are removed.
    #cdfpackagename: The cdf package name.  
    #probepackagename: The probe package name.
    
    require(cdfpackagename, character.only=TRUE)
    require(probepackagename, character.only=TRUE)
    
    probe.env.orig = get(probepackagename)
    
    #Remove probesets from the CDF environment
    rm(list=listOutProbeSets, envir=get(cdfpackagename))
    
    ##Set the PROBE env accordingly 
    ## (idea originally from gcrma compute.affinities.R)
    tmp = get('xy2indices', paste('package:', cdfpackagename, sep=''))
    
    newAB   = new('AffyBatch', cdfName=cleancdf)
    pmIndex =  unlist(indexProbes(newAB, 'pm'))
    subIndex = match(tmp(probe.env.orig$x, probe.env.orig$y, 
                         cdf=cdfpackagename), pmIndex)
    
    iNA = which(is.na(subIndex))
    
    ##Need to unlock the environment binding to alter the probes
    ipos = grep(probepackagename, search())
    env = as.environment(search()[ipos])
    
    unlockBinding(probepackagename, env) 
    assign(probepackagename, probe.env.orig[-iNA,], env=env)
    lockBinding(probepackagename, env )    
}


