library(optparse)
library(bnlearn)
library(quantmod)
library(infotheo)
library(igraph)


####################################################################################################
# Utility functions

option_list = list(
    make_option(c("-o", "--out"), type="character", default="multimodal_net/",
                help="output directory name [default= %default]", metavar="character"),
    make_option(c("-n", "--nodes"), type="integer", default=50,
                help="Number of nodes [default= %default]", metavar="integer"),
    make_option(c("-m", "--max_degree"), type="integer", default=3,
                help="Maximum degree of nodes [default= %default]", metavar="integer"),
    make_option(c("-d", "--discrete_prop"), type="double", default=0.2,
                help="Proportion of discrete nodes in the graph [default= %default]", metavar="double"),
    make_option(c("-s", "--n_samples"), type="integer", default=1000,
                help="Number of samples [default= %default]", metavar="integer"),
    make_option(c("-g", "--graph_file"), type="character", default=NULL,
                help="File of the underlying graph [default= %default]", metavar="character")
)

plotVariablesDistributionsToUniquePdf <- function(tbl, from, to, file){
    pdf(file)
    num = (to - from)/5
    if((to - from)%%5 != 0)
    num = num + 1
    par(mar = rep(2, 4), mfrow = c(num,5))

    ## plot freq variables
    for(i in from:to){
    name = colnames(tbl)[i]

    if(length(unique(tbl[,i])) < 30){
        counts <- table(tbl[,i])
        barplot(counts, main=name, xlab = "Values")
    } else {
        plot(density(tbl[,i]), main=name)
    }
    }
    dev.off()
}



####################################################################################################
# Node and interaction functions

gen_mixture = function(n, dist_list, flag=NULL){
    # Generates a multimodal distribution for n sample by drawing samples from a list of
    #distribution according to the flag vector which indicates from which distribution each sample
    #is to be drawn.
    #
    # Arguments :
    #    n            : Number of samples
    #    dist_list    : List of distributions of n samples
    #    flag        : Vector of length n with unique values 1:# of distributions
    #
    # Returns a multimodal continuous distribution.

    if(is.null(flag)){
        flag = sample(1:length(dist_list), size = n, replace = T) # Random probabilities
    }
    output_distribution = rep(0, length = n)

    i=1
    for(dist in dist_list){
        output_distribution = output_distribution + dist * (as.numeric(flag==i))
        i = i+1
    }
    return(output_distribution)
}

rescale = function(dist, min=0, max=1){
    # Rescale a given distribution to have values between min and max.
    #
    # Arguments :
    #    dist    : The continuous distribution to be rescaled.
    #    min,max : The range that the new distribution takes
    #
    # Returns the rescaled distribution.
    if(is.null(ncol(dist))){
        dist =    (dist-min(dist)) /(max(dist) - min(dist))
    }
    else{
        for(i in 1:ncol(dist)){
        dist[,i] = ( (dist[,i]-min(dist[,i])) /(max(dist[,i]) - min(dist[,i])) )
        }
    }
    dist = dist * (max-min) + min
    return(dist)
}

discretize_node = function(node_dist){
    # Discretize the distribution of a given node by either finding the valleys of its kernel
    # density estimation. If any single bin has more than 90% of the values and its standard
    # deviation is superior to a threshold (0.05), its content is re-discretized with an unsupervised
    # equal-frequencies discretization method into log(n_sample)-2 bins.
    # If the number of valleys is 0 or > 15, discretizes with the unsupervised method with
    # log(n_sample) bins.
    #
    # Arguments :
    #    node_dist    : The continuous distribution of a given node.
    #
    # Returns the discretized distribution.

    require(quantmod)
    require(infotheo)
    node_density = density(node_dist)
    # Use valleys in the density estimate as cut points.
    cut_points = node_density$x[quantmod::findValleys(node_density$y)]

    if((length(cut_points) != 0) & (length(cut_points<15))){
        # If we have found (not too many) valleys.
        node_levels = cut(node_dist, breaks = c(min(node_dist), cut_points, max(node_dist)),
                          labels = 1:(length(cut_points)+1), include.lowest = T)
        # Testing each resulting bin.
        uninformative_bins = sapply(unique(node_levels),
                                    function(x) (table(node_levels)[x] > (length(node_dist)*0.9)) && # bin has > 90% of values
                                                (sd(node_dist[node_levels==x]) > 0.05)) # and its content sd is > 0.05
        if(any(uninformative_bins)){
            uninformative_level = which(uninformative_bins)
            bin_levels = infotheo::discretize(node_dist[node_levels==uninformative_level],
                                                disc="equalwidth", nbins=log(length(node_dist))-2)[,1]

            # Merge new levels with previous node levels
            bin_levels = as.numeric(bin_levels) + uninformative_level - 1
            node_levels = as.numeric(node_levels)
            node_levels[node_levels > uninformative_level] = node_levels[node_levels > uninformative_level] +
                                                             max(bin_levels) - 1
            node_levels[node_levels == uninformative_level] = bin_levels
            node_levels = as.factor(node_levels)
        }
        node_levels = droplevels(node_levels) # Drop unused levels
        levels(node_levels) = 1:nlevels(node_levels) # Remap levels
    }
    else{ # If there are no valleys, perform equalwidth discretization
        node_levels = infotheo::discretize(node_dist, disc = "equalwidth", nbins = log(length(node_dist)))[,1]
    }
    if(median(table(node_levels)) < (length(node_dist)/100)){ # If peak and valleys result in too many almost empty bins, do equal width
        node_levels = infotheo::discretize(node_dist, disc = "equalwidth", nbins = log(length(node_dist)))[,1]
    }
    return(as.factor(node_levels))
}


gaussian_multivariate = function(parent_nodes_dist_df){
    # Generates a node distribution as the sum of its parents with gaussian noise.
    #
    # Arguments :
    #    node_dist    : The continuous distribution of a given node.
    #
    # Returns the child distribution.

    n_parents = dim(parent_nodes_dist_df)[2]
    n_sample = dim(parent_nodes_dist_df)[1]

    coeffs = runif(n_parents, min = 0.5, max=1)
    coeffs = t(matrix(rep(coeffs,n_sample), nrow=2))
    noise = matrix(rnorm(n_sample*n_parents), ncol=n_parents)

    child_dist = rowSums((coeffs*parent_nodes_dist_df) + 0.1*noise)
    return(rescale(child_dist))

}

linear_SEM = function(parent_nodes_dist_df){

    n_parents = dim(parent_nodes_dist_df)[2]
    n_sample = dim(parent_nodes_dist_df)[1]

    coeffs = sapply(1:n_parents, function(i) if(runif(1)>0.5){runif(1,-2,-0.1)} else{runif(1,0.1,2)})
    noise = matrix(rnorm(n_sample*n_parents, 0, 1), ncol=n_parents)
    K = matrix(rep(runif(n_parents, 0.1, 0.5), each=n_sample), ncol=n_parents)
    alpha = matrix(rep(runif(n_parents, 2, 4), each=n_sample), ncol=n_parents)
    noise = sign(noise) * abs(noise)**(alpha) * K

    child_dist = rowSums((coeffs*(parent_nodes_dist_df)) + noise)
    plot3d(X,Y,child_dist)

}

nonLinearSEM = function(node, parent_nodes_dist_df, discrete_parents){

    n_parents = dim(parent_nodes_dist_df)[2]
    n_sample = dim(parent_nodes_dist_df)[1]
    child_dist = numeric(n_sample)
    node_params[[node]] <<- list()
    continuous_parents = which(! 1:n_parents %in% discrete_parents)

    # build the full dataset with crossterms for computing resulting child distribution
    n_cross_terms = choose(length(continuous_parents), 2) # Add cross terms if enough parents, i.e. parent1*parent2
    coeffs = rep(1,n_parents+n_cross_terms)
    if(n_cross_terms>0){
        coeffs[1:n_parents] = 0 # Use only cross terms
        for(i in 1:(length(continuous_parents)-1)){
            for(j in (i+1):length(continuous_parents)){
                parent_nodes_dist_df = cbind(parent_nodes_dist_df, ((rescale(parent_nodes_dist_df[,continuous_parents[i]],-1,1)) *
                                                                    (rescale(parent_nodes_dist_df[,continuous_parents[j]],-1,1))))
            }
        }
    }

    polynomial_parents = setdiff(1:(n_parents+n_cross_terms), discrete_parents) #Continuous parents and crossterms have exponents
    for(parent in polynomial_parents){ 
        parent_nodes_dist_df[,parent] = rescale(parent_nodes_dist_df[,parent],-1,1) #Rescale between -1 and 1
    }

    if(length(discrete_parents)==0){ # Continuous parents -> continuous child
        # Draw random exponents
        exponents = sample(1:3, size = length(polynomial_parents), replace = T)
        for(parent in polynomial_parents){ 
            parent_nodes_dist_df[,parent] = parent_nodes_dist_df[,parent]**exponents[parent]
        }
        noise = matrix(rnorm(n_sample*(n_parents+n_cross_terms)), ncol=(n_parents+n_cross_terms))
        child_dist = rowSums( t(t(rescale(parent_nodes_dist_df)) * coeffs) + t( t(noise)*(0.1/n_parents/exponents)) )
    }

    else{ # Mixed parents -> continuous child
        # Draw different exponents for each unique observed combinations of discrete parents
        discrete_parent_levels = as.matrix(unique(parent_nodes_dist_df[,discrete_parents]),ncol=length(discrete_parents))
        comb_rows = numeric(nrow(parent_nodes_dist_df))
        for(comb_id in 1:nrow(discrete_parent_levels)){
            comb = discrete_parent_levels[comb_id, ]
            comb_rows[apply(as.matrix(parent_nodes_dist_df[,discrete_parents], ncol=length(discrete_parents)),
                            1, function(row) all(row == comb))] = comb_id
        }
        max_unique_combs = 3**(length(continuous_parents)) #Number of possible combinations according to the number of continuous parents
        if(nrow(discrete_parent_levels) > max_unique_combs){ #Merge some parent combinations together so that we have max_unique_combs unique
            reduced_comb_ids = sample(1:max_unique_combs, max_unique_combs, FALSE)
            reduced_comb_ids = c(reduced_comb_ids, sample(1:max_unique_combs, nrow(discrete_parent_levels)-max_unique_combs, TRUE))
            for(reduced_comb_id in 1:max_unique_combs){
                comb_rows = reduced_comb_ids[comb_rows]
            }
        }
        comb_ids = unique(comb_rows)
        set_of_exponents = expand.grid(rep(list(1:3), length(continuous_parents)))# Set of all possible exponents given the parents combinations
        node_params[[node]]$comb_rows <<- comb_rows
        node_params[[node]]$comb_ids <<- comb_ids
        node_params[[node]]$exponents_comb <<- list()
        node_params[[node]]$continuous_parents <<- continuous_parents

        for(i in 1:length(comb_ids)){
            comb_id = comb_ids[i]
            which_rows = (comb_rows==comb_id)
            # Draw new exponents
            #exponents = sample(1, size = length(continuous_parents)+n_cross_terms, replace = T)
            exponents = rep(0, size = n_parents+n_cross_terms)
            exponents[continuous_parents] = as.numeric(set_of_exponents[sample(nrow(set_of_exponents),1),])
            node_params[[node]]$exponents_comb[[comb_id]] <<- exponents[continuous_parents]

            parents_mix = matrix(0, nrow=length(which(which_rows)), ncol=n_parents+n_cross_terms)

            for(parent in polynomial_parents){
                #parents_mix[,parent] = (rescale(parent_nodes_dist_df[which_rows,parent],-1,1))**exponents[parent]
                parents_mix[,parent] = parent_nodes_dist_df[which_rows,parent]**exponents[parent]
                parents_mix[,parent] = parents_mix[,parent] + rnorm(length(which(which_rows)), sd=0.1)
            }
            child_dist[which_rows] = rowSums(as.matrix(rescale(parents_mix[,continuous_parents])))
            #child_dist[which_rows] = rowSums(as.matrix(parents_mix[,continuous_parents]))
            set_of_exponents = set_of_exponents[apply(set_of_exponents, 1, function(row) !all(row==exponents[continuous_parents])),,drop=F]
        }
    }
    #plot3d(X,Y,child_dist, col = colors)
    return(rescale(child_dist))
    #return(rescale(cos(rescale(child_dist)*5)))
    #return(rescale((rescale(child_dist)*2-1)**2))
}


syntren_like = function(node, parent_nodes_dist_df){
    # Compute the distribution of a node given its parents' distributions with Michaelis-Menten
    #and Hill kinetic models. This function draws random reaction parameters and calls mm_kinetics().
    # This implementation is based on the equations used in [1] and [2].
    # [1] Van den Bulcke, Tim, et al. "SynTReN: a generator of synthetic gene expression data for design and analysis of structure learning algorithms." BMC bioinformatics 7.1 (2006): 43.
    # [2] Mendes, Pedro, Wei Sha, and Keying Ye. "Artificial gene networks for objective comparison of analysis algorithms." Bioinformatics 19.suppl_2 (2003): ii122-ii129.
    #
    # Arguments :
    #    parent_nodes_dist_df    : A dataframe that contains the distributions of the node's parents.
    #
    # Returns the node's distribution.

    # Assign parents either an inhibitor or activator role.
    n_parents = dim(parent_nodes_dist_df)[2]
    n_sample = dim(parent_nodes_dist_df)[1]
    P = sample(x = 1:n_parents, size = 1) # number of activators
    Q = n_parents - P # number of inhibitors

    A_columns = sample(1:n_parents, P, replace = F)

    # Retrieve the Activators and Inhibitors concentrations
    A = as.matrix(parent_nodes_dist_df[, A_columns])
    if(Q > 0){
        I_columns = setdiff(1:n_parents, A_columns)
        I = as.matrix(parent_nodes_dist_df[, I_columns])
    }
    else{
        Q = 1
        I = matrix(rep(0, n_sample), ncol = 1) # Add an inhibitor with [I] = 0
    }

    # Draw a basal rate of reaction.
    V0 = runif(1, 0, 0.5) # Basal output
    Vi = runif(P, 0.5, 1) # Activators' Vmax

    # The variables K represent concentrations at which the effect of the inhibitor or activator
    #is half of its saturating value.
    Ki = rnorm(P, 0.5, 0.15)
    Kj = rnorm(Q, 0.5, 0.15)

    # The hill coefficient n describes the fraction of the macromolecule saturated by ligand as a
    #function of the ligand concentration; it is used in determining the degree of cooperativeness.
    #    - n>1 - Positively cooperative binding: Once one ligand molecule is bound to the enzyme,
    #        its affinity for other ligand molecules increases.
    # In the Syntren paper, example interactions are given for n = 1,2,5,10.
    # If there is one parent, n is either 1 or 2 so that the interaction function is not too steep.
    # If there are more parents, ns are drawn such that abs(n_k-n_l) < 2 for all k,l with 1 < n < 10.
    # That way, no single interaction is masked by a much stronger interaction with one of the other
    # parents and the sampled distributions are faithful to the graph. 
    if(n_parents > 2){
        N = rbinom(1, 10, 0.2)+1
    }
    else{
        N = 1
    }
    ni = round(runif(P, max(N-1.49, 0.51), N+1.49))
    nj = round(runif(Q, max(N-1.49, 0.51), N+1.49))

    node_params[[node]] <<- list(A=A, P=P, I=I, Q=Q, V0=V0, Vi=Vi, ni=ni, nj=nj, Ki=Ki, Kj=Kj, A_columns=A_columns)

    noise = scale(rgamma(n_sample, shape=2, rate=20), scale = F) * sample(c(-1,1), size = n_sample, replace = T)
    child_dist = rescale(mm_kinetics(A, P, I, Q, V0, Vi, ni, nj, Ki, Kj)) + noise
    child_dist = apply(child_dist, 1, function(x) max(0,x))

    return(rescale(child_dist))
}


mm_kinetics = function(A, P, I, Q, V0 = NULL, Vi=NULL, ni = NULL, nj = NULL, Ki = NULL, Kj = NULL){
    # Compute a distribution with the general steady-state equation for P activator and Q inhibitors
    #based on Michaelis-Menten and Hill kinetics.
    #
    # Arguments :
    #    A    : Matrix with the activators concentrations (dim : nsample * P).
    #    P    : Number of activators.
    #    I    : Matrix with the inhibitors concentrations (dim : nsample * Q).
    #    Q    : Number of inhibitors.
    #    V0    : Basal rate of the output production.
    #    Vi    : Vector with the Vmax of the activators.
    #    ni    : Vector with the Hill constants of the activators.
    #    nj    : Vector with the Hill constants of the inhibitors.
    #    Ki    : Vector with the Michaelis constants of the activators.
    #    Kj    : Vector with the Michaelis constants of the inhibitors.
    #
    # Returns the output distribution.

    dividend = V0
    for(i in 1:P){
        product = 1
        for(j in 1:P){
            if(j!=i){
                product = product * (1 + (A[,j]/Ki[j])^ni[i] )
            }
        }
        product = product * Vi[i]
        dividend = dividend + (A[,i]/Ki[i])^ni[i] * product
    }
    product = 1
    for(i in 1:P){
        product = product * (1 + (A[,i]/Ki[i])^ni[i])
    }
    divisor = product
    product = 1
    for(j in 1:Q){
        product = product * (1 + (I[,j]/Kj[j])^nj[j])
    }
    divisor = divisor * product

    output_distribution = dividend / divisor
    return(output_distribution)
}



#draw_child_dist = function(parent_nodes_dist_df){
#    # Draws a random sample from parents' distributions with the same scheme as the multimodal
#    #distribution generation method.
#
#    if(!is.null(dim(parent_nodes_dist_df))){ # if more than 1 parent
#        n_parents = dim(parent_nodes_dist_df)[2]
#        n_sample = dim(parent_nodes_dist_df)[1]
#        probs = runif(n = n_parents)
#        probs = probs/sum(probs)
#        flag = sample(1:n_parents, size = n_sample, replace = T, prob = probs)
#
#        #flag = sapply(X=unique(flag), FUN = function(x) as.numeric(flag == x))
#        #inverse_flag = abs(1-flag)
#        #participative_noise = runif(n = n_sample, min = 0, max = 0.5)
#        #m = (matrix(runif(n_sample * n_parents), nrow = n_sample) * inverse_flag)
#        #m = m/apply(m, 1, sum) * participative_noise
#        #participative_noise = flag*(1-participative_noise)
#        #weighted_flag = participative_noise + m
#
#        X = rep(0, n_sample)
#        for(i in 1:n_parents){
#            X = X + parent_nodes_dist_df[,i] * as.numeric(flag == i) + runif(n_sample, -0.1, 0.1)
#            #X = X + parent_nodes_dist_df[,i] * weighted_flag[,i] + runif(n_sample, -1, 1) * weighted_flag[,i]
#        }
#    }
#    else{ # if only one parent
#        noise = scale(rgamma(length(parent_nodes_dist_df), shape=2, rate=50), scale = F) *
#                sample(c(-1,1), size = length(parent_nodes_dist_df), replace = T)
#        X = parent_nodes_dist_df + noise #runif(n = length(parent_nodes_dist_df), min = -0.1, max=0.1)
#    }
#    return(X)
#}

distFromDiscrete <- function(parent_nodes_dist_df, nBins){
    probTable = data.frame(parent_nodes_dist_df[!duplicated(parent_nodes_dist_df),])
    probTableList=list()
    for(i in 1:nrow(probTable)){
    weights <- runif(nBins)
    f <- 1/sum(weights)
    weights <- f*weights
    weights <- weights * runif(nBins,0,2)
    weights <- weights/sum(weights)
    probTableList[[paste(probTable[i,], collapse = "-")]] = cumsum(weights)
    }

    vec = rep(NA, nrow(parent_nodes_dist_df))
    for(i in 1:nrow(parent_nodes_dist_df)){
    r = runif(1)
    vecProb = probTableList[[paste(parent_nodes_dist_df[i,], collapse = "-")]]
    for(col in 1:length(vecProb)){
        if(r <= vecProb[col]){
        vec[i] = col
        break
        }
    }
    }
    return(as.factor(vec))
}


#continualize_node = function(node_dist){
#
#    number_of_modes = nlevels(node_dist)
#
#    peaks = seq(-10, 10, by=min(5, 20/number_of_modes))
#    peaks = peaks + rnorm(length(peaks), mean = 0, sd=0.5)
#    means = sample(peaks, size = number_of_modes)
#    means = sort(means)
#    sds = runif(number_of_modes, min = 1, max = 2)
#
#    dist_list = lapply(X = 1:number_of_modes,
#                        FUN = function(x) rnorm(n=length(node_dist), mean = means[x], sd = sds[x]))
#    node_dist_continuous = rep(0, length(node_dist))
#
#    for(level in levels(node_dist)){
#        n_sample_level = length(which(node_dist==level))
#        probs = runif(number_of_modes, 0, min(0.2, 1/number_of_modes))
#        level_index = which(levels(node_dist) == level)
#        probs[level_index] = 1 - sum(probs) + probs[level_index]
#        flag = sample(1:number_of_modes, size=n_sample_level, replace = T, prob = probs)
#        node_dist_continuous[node_dist==level] = gen_mixture(n_sample_level, dist_list, flag)
#    }
#    return(rescale(node_dist_continuous))
#}


draw_values_for_discrete = function(node_dist){
    # Draw real values between 0 and 1 for eachh level of a discrete distribution. The information
    # between levels is randomized, i.e. levels "1" and "2" are not guaranteed to have real values
    # that are closer than those of "1" and "3".
    #
    # Arguments :
    #    node_dist    : A vector of factors
    #
    # Returns a vector of numerical values between 0 and 1, where length(unique(output)) = nlevels(input)

    node_levels = levels(node_dist)
    step_size = 1/(nlevels(node_dist)+1)
    real_values = seq(step_size, 1-step_size, by = step_size) + rnorm(nlevels(node_dist), 0, step_size/5)
    real_values = sample(real_values) # Randomize
    node_dist = as.numeric(node_dist)
    for(i in 1:length(node_levels)){
        node_dist[node_dist == node_levels[i]] = real_values[i]
    }
    return(node_dist)
}


generate_child_dist = function(node, parents, parent_nodes_dist_df){
    # Generates a node's distribution from its parents' by calling various subfunctions.
    # Here is the general scheme :
    # - If all of the node's parents are discrete, it becomes discrete itself and its distribution
    #    is drawn with the fonction distFromDiscrete, with a number of bins between 2 and 5.
    # - else, if the node is discrete : discretize its parents if needed and draw values with distFromDiscrete for 
    #                                   a number of bins between 2 and 5.
    #         if the node is continuous : run non-linear structural equation models with crossterms if all parents
    #                                     are continuous, or in polynomial form where exponents depend on the state
    #                                     of the discrete parents if there are any, with function nonLinearSEM.
    #
    # Arguments :
    #    node                    : A string which is the label of the node.
    #    parents                 : A vector of strings for the node's parents.
    #    parent_nodes_dist_df    : A dataframe containing the parent nodes distribution.
    #
    # Returns the node distribution

    if(all(as.logical(is_node_discrete[parents]))){
        nBins = 2#sample(2:5,1)
        node_dist = distFromDiscrete(data.frame(parent_nodes_dist_df), nBins=nBins)
        node_params[[node]] <<- list(discrete=TRUE, parents=parents, method="distFromDiscrete", nBins=nBins)
        is_node_discrete[node] <<- TRUE
    }
    else{
        discrete_parents = parents[as.logical(is_node_discrete[parents])]
        if(is_node_discrete[[node]]){
            for(continuous_parent in setdiff(parents, discrete_parents)){
                parent_nodes_dist_df[[continuous_parent]] = discretize_node(parent_nodes_dist_df[[continuous_parent]])
                #parent_nodes_dist_df[[continuous_parent]] = infotheo::discretize(parent_nodes_dist_df[[continuous_parent]], disc = "equalwidth", nbins = 3)[,1]
            }
            nBins = 2#sample(2:5,1)
            node_dist = distFromDiscrete(data.frame(parent_nodes_dist_df), nBins=nBins)
            node_params[[node]] <<- list(discrete=TRUE, parents=parents, method="distFromDiscrete", method_parents="discretized", nBins=nBins)
        }
        else{
            node_dist = nonLinearSEM(node, data.frame(parent_nodes_dist_df), which(as.logical(is_node_discrete[parents])))
            node_params[[node]] <<- c(node_params[[node]], discrete=FALSE, parents=parents, method="nonLinearSEM", discrete_parents=is.null(discrete_parents))
        }
    }
    return(node_dist)

}



####################################################################################################
# Run the simulation


# Retrieve command line args
args = commandArgs(trailingOnly=TRUE)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Generate network
nodes = paste0("X", 1:opt$nodes)
g = random.graph(nodes, method="melancon", max.degree=opt$max_degree) # Draws a random DAG with the uniform method.
#g = random.graph(nodes, method="ordered", prob = 4/length(nodes))
#if(is.null(opt$graph_file)){
#    g = random.graph(nodes, method="melancon", max.degree=opt$max_degree)
#} else{
#    g = read.net(opt$graph_file)
#}
n = opt$n_samples # Number of samples.
simulation = as.data.frame(x = matrix(nrow = n, ncol = length(nodes)))
colnames(simulation) = nodes

discrete_prob = opt$discrete_prop # The approximative percentage of discrete nodes.
discrete_flag = discrete_prob>=(1/n)
# Handle discretization
is_node_discrete = setNames(as.list(rep(FALSE, length(nodes))), nodes) #List/dictionary of nodes.
if(discrete_flag){
    #nodes_to_discretize = nodes[sample(c(TRUE,FALSE), length(nodes), TRUE, c(discrete_prob, 1-discrete_prob))]
    nodes_to_discretize = sample(nodes, floor(discrete_prob*length(nodes)), FALSE)
    is_node_discrete[nodes_to_discretize] = TRUE
}


orphan_nodes = nodes[!nodes %in% g$arcs[,2]] # Nodes without parents.
for(node in orphan_nodes){
    # Generate a multimodal distribution for each orphan node
    number_of_modes = floor(runif(1, min = 2, max = 5)) # modes in the node's distribution
    peaks = seq(-20, 20, by=5) # 9 peaks to choose from
    peaks = peaks + rnorm(length(peaks), mean = 0, sd=0.5) # Add randomness to the choice of peaks
    means = sample(peaks, size = number_of_modes) # means of the individual distributions (modes of the multimodal distribution)
    sds = runif(number_of_modes, min = 1, max = 1) # standard errors of the individual distributions

    dist_list = lapply(X = 1:number_of_modes, FUN = function(x) rnorm(n, mean = means[x], sd = sds[x]))

    if(is_node_discrete[[node]]){
        #simulation[,node] = sample(1:sample(2:2,1), size = n, replace = T)#discretize_node(simulation[,node])
        simulation[,node] = sample(1:2, size = n, replace = T)#discretize_node(simulation[,node])
    }
    else{
        simulation[,node] = gen_mixture(n, dist_list) # generates the multimodal distribution from the constructed normal distributions
        simulation[,node] = rescale(simulation[,node]) # rescale to [0,1]
    }
}


# Simulate the rest of the network, following the causality order.
node_params = list() # Where the simulation parameters are saved for each node.
simulation_finished = F
list_of_sampled_nodes = orphan_nodes # List of nodes which currently have a distribution.

while(!simulation_finished) {
    # For every node which has no distribution yet
    for(node in nodes[!nodes %in% list_of_sampled_nodes]){
        # Find its parents
        node_parents = g$arcs[,1][g$arcs[,2] == node]
        # If all of its parents have been sampled
        if(all(node_parents %in% list_of_sampled_nodes)){
            dist_parents = simulation[node_parents]
            simulation[,node] = generate_child_dist(node, node_parents, dist_parents) #Generate distribution of given node
            list_of_sampled_nodes = c(list_of_sampled_nodes, node)
        }
    }
    simulation_finished = all(nodes %in% list_of_sampled_nodes)
}


####################################################################################################
# Save and plot

dir.create(file.path(opt$out), showWarnings = FALSE)

save(node_params, file=paste0(opt$out, "/node_params.RData"))
write.table(simulation, file = paste0(opt$out, "/simulation.tsv"), row.names = F, quote = F, sep = '\t')
write.table(g$arcs, file = paste0(opt$out, "/true_edges.txt"), row.names = F, col.names = F, sep = '\t', quote = F)
g_ig = igraph::read.graph(format = "ncol", file = paste0(opt$out, "/true_edges.txt"), directed=T)
pdf(file=paste0(opt$out, "/graph.pdf"))
plot(g_ig, vertex.color = ifelse(is_node_discrete[V(g_ig)$name], 'skyblue', 'orange'),
     layout = layout_(g_ig, with_gem()))

for(i in 1:(((length(nodes)-1)%/%50)+1) ){
    start = ((i-1)*50)+1
    stop = min(length(nodes),i*50)
    plotVariablesDistributionsToUniquePdf(simulation, start, stop,
                                            paste0(opt$out, "/sim_dists_N-", start, "-", stop, ".pdf"))
}
