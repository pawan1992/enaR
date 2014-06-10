enaUncertainty=function(x = 'network object', lower.percent=0.9, upper.percent=1.1, iter=10){
	
   if (class(x) != 'network'){warning('x is not a network class object')}					   # check object class
	
	U = unpack(x)                                                      		       	        # unpack network object
	fluxes = which(U$F!=0, arr.ind=TRUE)                               		      	        # identify internal fluxes (from,to)
	inputs = seq(from=1, to=length(U$z), by=1)                         		        	        # identify inputs
	outputs = seq(from=1, to=length(U$y), by=1)                        		        		    # identify outputs
	vertex.names = rep(0, x$gal$n)													       # get names
	respiration=rep(0, x$gal$n)
	living=U$Living
	storage=U$X
	for(i in 1:x$gal$n){
		vertex.names[i]=x$val[[i]]$vertex.names
		}
	
	E = matrix(0, nrow=x$gal$n, ncol=(nrow(fluxes)+length(inputs)+length(outputs))) 		    # initialize E
	
	for(i in 1:length(U$z)){                                                         	         # add inputs to E
		E[i,i] = U$z[i]
		}

	for(i in 1:nrow(E)){																# add outputs to E
		E[i,length(U$z)+outputs[i]] = -U$y[i]
		}
	
	for(f in 1:nrow(fluxes)){															# add internal node inputs
		E[fluxes[f,2],(length(U$z)+length(U$y)+f)] = U$F[fluxes[f,1],fluxes[f,2]]
		}
		
	for(f in 1:nrow(fluxes)){															# add internal node outputs
		E[fluxes[f,1],(length(U$z)+length(U$y)+f)] = -U$F[fluxes[f,1],fluxes[f,2]]
		}
		
	F = rep(0,x$gal$n)																	# create F
	
	G = rbind(diag(rep(1, ncol(E))), diag(rep(-1, ncol(E))))								# create G
	
	if(length(lower.percent)==1 & length(upper.percent)==1){
		
	H = c(rep(lower.percent, ncol(E)),rep(-upper.percent, ncol(E)))							# create H
	}else{
	H=c(lower.percent,-upper.percent)
	}
	
	xs = xsample(E=E, F=F, G=G, H=H, iter=iter)										# calculate plausible coefficients
	
	
	z.ena = rep(0, x$gal$n)															# initialize
	y.ena = rep(0, x$gal$n)
	F.ena = matrix(0, nrow=nrow(U$F), ncol=ncol(U$F))
	plausible.models = list()
	
	for(k in 1:nrow(xs$X)){															# construct pluasible models
		for(z in 1:length(U$z)){														# inputs (z)
			z.ena[z] = xs$X[k,z]*U$z[z]
		}
		for(y in (length(U$z)+1):(length(U$z)+length(U$y))){								# outputs (y)
			y.ena[y-length(U$z)] = xs$X[k,(y)]*U$y[y-length(U$z)]
		}
		for(f in (length(U$z)+length(U$y)+1):ncol(xs$X)){									# internal fluxse (F)
			F.ena[fluxes[(f-length(U$z)-length(U$y)),1], fluxes[(f-length(U$z)-length(U$y)),2]] = xs$X[k,f]*U$F[fluxes[(f-length(U$z)-length(U$y)),1],fluxes[(f-length(U$z)-length(U$y)),2]]
		}
		
		rownames(F.ena) = vertex.names
		colnames(F.ena) = vertex.names
		
		plausible.models[[k]] = pack(flow=F.ena, input=z.ena, export=y.ena, respiration=respiration, living=living, storage=storage)
	}
	
	return(plausible.models)
}
