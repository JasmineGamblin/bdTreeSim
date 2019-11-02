## State Update for Gillespie Algorithm

gillespie_update <- function(lambda, mu, t, N) {
  #compute rates
  T_birth <- lambda*N
  T_death <- mu*N

  #update t
  t <- t + rexp(1, rate = T_birth + T_death)

  u <- runif(1, min = 0, max = T_birth + T_death)

  #birth event
  if (u < T_birth) {
    N <- N+1
  }

  #death event
  else {
    N <- N-1
  }
  return(list("t" = t, "N" = N))
}



## compute the value N_s of the trajectory (t,N) at time s

traj_N <- function(traj, s) {
  idx <- max(which(traj$t<=s))
  return(traj$N[idx])
}

#' @title Trajectory Simulation
#'
#' @description Simulates a trajectory of the population size under a constant-rate birth-death model as well as the removal vector \eqn{y} (\eqn{y_i=1} iff \eqn{s_i} is removed upon sampling). It uses Gillespie algorithm to simulate portions of trajectory in the time intervals between removals.
#'
#' @param lambda birth rate
#'
#' @param mu death rate
#'
#' @param r removal probability, the probability for an individual to be removed upon sampling
#'
#' @param sampling_times vector containing the sampling times \eqn{(s_j)_{j=1...k}} in increasing order
#'
#' @return an object of the class \code{trajectory}, containing:
#' \itemize{
#'    \item \code{t} a sequence of time points \eqn{(t_i)_{i=1...n}} with \eqn{t_1=0} and \eqn{t_n=s_k},
#'    \item \code{N} a sequence of population size values where \eqn{N_i} is the value on the interval \eqn{[t_i,t_{i+1})}, and
#'    \item \code{removed} a vector with the sampling times of removed samples.
#' }
#'
#' @examples trajectory(lambda=2, mu=1, r=0.2, sampling_times=c(1,2))
#'
#' @importFrom stats runif rexp
#'
#' @export

trajectory <- function(lambda, mu, r, sampling_times) {
  n_s <- length(sampling_times)
  removed <- sampling_times[runif(n_s, 0, 1) < r]
  N <- 1
  t <- 0
  i <- 1

  for (s in unique(c(removed, sampling_times[n_s]))) {
    while (t[i] < s && N[i] > 0) {
      update <- gillespie_update(lambda, mu, t[i], N[i])
      t <- c(t, update$t)
      N <- c(N, update$N)
      i <- i+1
    }
    if (N[i] == 0 && t[i] < s) {
      break
    }
    else {
      t[i] <- s
      N[i] <- N[i-1] - if (s %in% removed) 1 else 0
    }
  }

  traj <- list("t" = t, "N" = N, "removed" = removed)
  class(traj) <- 'trajectory'
  return(traj)
}


#' @title Weight of a Trajectory
#'
#' @description Computes the log probability density of the sampling times given the trajectory.
#'
#' @param traj object of class \code{trajectory}
#'
#' @param sampling_times vector containing the sampling times \eqn{(s_j)_{j=1...k}} in increasing order
#'
#' @param psi sampling rate
#'
#' @return a number
#'
#' @examples weight(traj=list('t'=c(0,0.5,1,2), 'N'=c(1,2,1,1), 'removed'=c(1)),
#' sampling_times=c(1,2), psi=0.3)
#'
#' @export

weight <- function(traj, sampling_times, psi) {

  if (traj$t[length(traj$t)] < sampling_times[length(sampling_times)]) {
    return(-Inf)
  }
  else {
    interval_widths <- diff(traj$t)
    interval_Ns <- traj$N[-length(traj$N)]
    logp <- -sum(interval_widths*interval_Ns)*psi

    for (s in sampling_times) {
      if (s %in% traj$removed) {
        logp <- logp + log((traj_N(traj, s) + 1)*psi)
      }
      else {
        logp <- logp + log(traj_N(traj, s)*psi)
      }
    }
    return(logp)
  }
}



## expand the length of the root branch of a tree

expand <- function(tree, t) {
  tree$branch_length <- tree$branch_length + t
  return(tree)
}


#' @title Tree Simulation from a Trajectory
#'
#' @description Reconstructs a tree from a trajectory and known sampling times.
#'
#' @param traj an object of class \code{trajectory}
#'
#' @param sampling_times vector containing the sampling times \eqn{(s_j)_{j=1...k}} in increasing order
#'
#' @return an object of class \code{tree}, containing an incident branch length \code{branch_length} and a list of 0 to 2 subtrees \code{subtrees}
#'
#' @examples tree_from_trajectory(traj=list('t'=c(0,0.5,1,2), 'N'=c(1,2,1,1),
#' 'removed'=c(1)), sampling_times=c(1,2))
#'
#' @importFrom stats runif
#'
#' @export

tree_from_trajectory <- function(traj, sampling_times) {
  t <- sort(unique(c(traj$t, sampling_times))) # all events (start, birth and death, sampling)
  n_t <- length(t)

  n_lineages <- 1 # n_lineages = length(lineages)
  lineages <- list(list(branch_length = t[n_t]-t[n_t-1], subtrees = NULL))

  for (i in (n_t-1):2) {
    N_i <- traj_N(traj, t[i])
    N_i0 <- traj_N(traj, t[i-1])

    # sampling event with removal
    if (t[i] %in% traj$removed) {
      # create new lineage
      lineages <- c(lapply(lineages, function(x) expand(x, t[i]-t[i-1])), list(list(branch_length = t[i]-t[i-1], subtrees = NULL)))
      n_lineages <- n_lineages + 1
    }

    # sampling event without removal
    else if (t[i] %in% sampling_times) {
      if (runif(1) < n_lineages/N_i) {
        # select lineage to add sampling event
        k <- sample(1:n_lineages, 1)
        lineages <- c(lapply(lineages[-k], function(x) expand(x, t[i]-t[i-1])), list(list(branch_length = t[i]-t[i-1], subtrees = lineages[k])))
      }
      else {
        # create new lineage
        lineages <- c(lapply(lineages, function(x) expand(x, t[i]-t[i-1])), list(list(branch_length = t[i]-t[i-1], subtrees = NULL)))
        n_lineages <- n_lineages + 1
      }
    }

    # birth event
    else if (N_i > N_i0 && runif(1) < choose(n_lineages,2)/choose(N_i,2)) {
      # coalesce two lineages
      pair <- sample(1:n_lineages, 2, replace = F)
      lineages <- c(lapply(lineages[-pair], function(x) expand(x, t[i]-t[i-1])), list(list(branch_length = t[i]-t[i-1], subtrees = lineages[pair])))
      n_lineages <- n_lineages - 1
    }

    # expand all lineages
    else {
      lineages <- lapply(lineages, function(x) expand(x, t[i]-t[i-1]))
    }
  }

  tree <- lineages[[1]]
  class(tree) <- 'tree'
  return(tree)
}



#' @title Estimated Sample Size
#'
#' @description Returns the estimated number of independent samples based on the weights, using Rényi's entropy: \eqn{\frac{1}{\sum_i w_i^2}}.
#'
#' @param weights a vector of weight values associated to the sampled trajectories
#'
#' @param log boolean, \code{TRUE} if log weights are provided
#'
#' @return a number
#'
#' @examples ess_value(weights=c(3,0,456,1,0,0,5,57), log=FALSE)
#'
#' @export

ess_value <- function(weights, log) {
  if (all(is.infinite(weights))) {
    return(0)
  }
  else if (log) {
    exp_weights <- exp(weights)
    exp_weights <- exp_weights/sum(exp_weights)
    return(1/sum(exp_weights**2))
  }
  else {
    return(1/sum(weights**2))
  }
}



#' @title Tree Sampling
#'
#' @description Samples trees under a birth death model conditioned on sampling times.
#' \itemize{
#'     \item If \code{k1} is specified then the function simulates \code{k1} weighted trajectories, then subsample \code{k2} trajectories from them.
#'     \item If \code{e} is specified the function will simulate just enough weighted trajectories to obtain an ESS value greater or equal to \code{e}, and then subsample \code{k2} trajectories from them.
#' }
#' In both cases, the function then reconstructs one tree per trajectory and returns these trees.
#' One of \code{k1} or \code{e} has to be specified. If both are specified, \code{e} will be ignored.
#'
#' @param k1 number of trajectories to simulate, from which \code{k2} trajectories will be subsampled
#'
#' @param e ESS value to achieve
#'
#' @param k2 number of trees to simulate
#'
#' @param sampling_times vector containing the sampling times \eqn{(s_j)_{j=1...k}} in increasing order
#'
#' @param lambda birth rate
#'
#' @param mu death rate
#'
#' @param psi sampling rate
#'
#' @param r removal probability
#'
#' @return a list of \code{tree} objects
#'
#' @examples sample_trees(k1=NULL, e=100, k2=100, sampling_times=c(1,2), lambda=2, mu=1, psi=0.3, r=0.2)
#'
#' @export

sample_trees <- function(k1=NULL, e=NULL, k2, sampling_times, lambda, mu, psi, r) {
  if (!is.null(k1)) {
    # simulate k1 weighted trajectories
    weighted_trajs <- lapply(1:k1, function(x) trajectory(lambda, mu, r, sampling_times))
    weights <- sapply(weighted_trajs, function(x) weight(x, sampling_times, psi))
    print(c('zero weights proportion', length(which(weights == -Inf))/k1))

    # reduce size of weights to avoid overflow
    weights <- weights - max(weights)

    # subsample k2 trajectories from them
    exp_weights <- exp(weights)
    S <- sum(exp_weights)
    trajs <- sample(1:k1, k2, replace = T, prob = exp_weights)
    print(c('ESS value', ess_value(exp_weights/S, log=F)))
  }

  else if (!is.null(e)) {
    ess <- 0
    weighted_trajs <- c()
    weights <- c()

    # sample trajectories until the wanted ESS value is achieved
    while (ess < e) {
      weighted_trajs <- c(weighted_trajs, list(trajectory(lambda, mu, r, sampling_times)))
      weights <- c(weights, weight(weighted_trajs[[length(weighted_trajs)]], sampling_times, psi))
      ess <- ess_value(weights, log=T)
    }

    # reduce weights
    weights <- weights - max(weights)

    # subsample k2 trajectories
    trajs <- sample(1:length(weighted_trajs), k2, replace = T, prob = exp(weights))
    print(c('simulated trajectories', length(weighted_trajs)))
    print(c('ESS value', round(ess, digits =0)))
  }

  # simulate trees from the k2 trajectories
  trees <- lapply(weighted_trajs[trajs], function(x) tree_from_trajectory(x, sampling_times))
  return(trees)
}



#' @title Tree Height
#'
#' @description Compute the height of a tree, that is the time elapsed between the most recent common ancestor of the samples and the last sample.
#'
#' @param tree an object of class \code{tree}
#'
#' @param t_max the last sampling time
#'
#' @return a number
#'
#' @examples tree_height(tree=list('branch_length'=1, 'subtrees'=list(list('branch_length'=1,
#' 'subtrees'=NULL))), t_max=2)
#'
#' @export

tree_height <- function(tree, t_max) {
  return(t_max - tree$branch_length)
}




#' @title Compute the Height Distribution
#'
#' @description Compute the height distribution for a list of trees. The height of a tree is the time elapsed between the most recent common ancestor of the samples and the last sample.
#'
#' @param trees a list of \code{tree} objects
#'
#' @param t_max time of the last sampling event
#'
#' @param print boolean, if \code{TRUE} prints the median, mean and 95\% HPD interval for the distribution
#'
#' @param plot boolean, if \code{TRUE} plot the histogram of the distribution
#'
#' @param save boolean, if \code{TRUE} save the height distribution in a .RData file
#'
#' @param filename if \code{save} is \code{TRUE}, name of the file who will contain the distribution
#'
#' @return a vector containing the heights of the trees specified in \code{trees}
#'
#' @examples height_distribution(trees=list(list('branch_length'=1,
#' subtrees=list(list('branch_length'=1, subtrees=NULL)))), t_max=2, print=TRUE, plot=TRUE)
#'
#' @importFrom stats median
#'
#' @importFrom graphics hist
#'
#' @importFrom HDInterval hdi
#'
#' @export

height_distribution <- function(trees, t_max, print = FALSE, plot = FALSE, save = FALSE, filename = "height_dist.RData") {
  dist <- sapply(trees, function(x) tree_height(x, t_max))
  if (print) {
    print(c("median:", median(dist)))
    print(c("mean:", mean(dist)))
    print(c("95% HPD interval:", hdi(dist)))
  }
  if (plot) {
    hist(dist)
  }
  if (save) {
    save(dist, file = filename)
  }
  return(dist)
}



#' @title Write Tree in Newick Format
#'
#' @description Write the Newick description of a tree, and return the corresponding string.
#'
#' @param tree an object of type \code{tree}
#'
#' @return a string
#'
#' @examples str_tree(tree=list('branch_length'=1, 'subtrees'=list(list('branch_length'=1, 'subtrees'=NULL))))
#'
#' @export

str_tree <- function(tree) {
  if (is.null(tree$subtrees)) {
    return(paste(":",toString(tree$branch_length), sep = ""))
  }
  else if (length(tree$subtrees) == 1) {
    return(paste("(", str_tree(tree$subtrees[[1]]), "):", toString(tree$branch_length), sep = ""))
  }
  else {
    return(paste("(", str_tree(tree$subtrees[[1]]), ",", str_tree(tree$subtrees[[2]]),"):", toString(tree$branch_length), sep = ""))
  }
}



#' @title Save Trees in Newick Format
#'
#' @description Write a list of trees in Newick format in the specified text file.
#'
#' @param trees list of \code{tree} objects
#'
#' @param filename name of the file to save the trees
#'
#' @examples newick(trees=list(list('branch_length'=1, subtrees=list(list('branch_length'=1, subtrees=NULL)))),
#' filename='one_tree.txt')
#'
#' @export

newick <- function(trees, filename) {
  str_trees <- sapply(trees, str_tree)
  writeLines(str_trees, con = filename, sep = ";\n")
}








