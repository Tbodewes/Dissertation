library(bnlearn)
library(parallel)
library(stringr)
library(prodlim)
library(dplyr)
source("Scoring.R")


#' Initialize dataframe of possible moves for local search
#'
#' The local search algorithm represents possible moves as ordered pairs of
#' distinct nodes. Whether such a pair represents an arc addition, deletion or
#' reversal depends on whether the arc is absent (addition), present (deletion)
#' or reversed (reverse) in the dag that the move acts on.
#'
#' Moves are stored in a dataframe with as columns the from-node, the to-node, a
#' logical for whether the score difference for this move is still correct for
#' the current DAG, and the score difference for the move
#'
#' @param node.names character vectors of node names
#' @param blacklist dataframe indicating which directed arcs are not allowed in
#'   the DAG that is being learned. Consists of a 'from' and 'to' column
#' @param whitelist dataframe indicating which directed arcs must be in the DAG
#'
#' @return dataframe containing each move, its score and whether the score is
#'   up-to-date
#' @export
#'
#' @examples
initializeMoves <- function(node.names, blacklist, whitelist){
  #Get all pairs of distinct nodes
  moves <- expand.grid(from = node.names, to = node.names, stringsAsFactors = F)
  moves <- with(moves, moves[from != to,])
  
  #Efficiently enforce black- and whitelist by removing moves
  require(prodlim)
  if(!is.null(blacklist)){
    #Prevent arc additions
    names(blacklist) <- c("from", "to")
    moves <- moves[-na.omit(row.match(data.frame(blacklist), moves)),]
  }
  if(!is.null(whitelist)){
    #Prevent arc deletions and reversals
    whitelistReverse <- whitelist[, c(2,1)]
    names(whitelist) <- c("from", "to")
    names(whitelistReverse) <- c("from", "to")
    moves <- moves[-na.omit(row.match(rbind(whitelist, whitelistReverse), moves)),]
  }
  
  #Indicates that none of the scores for the moves has been computed, hence
  #all should be evaluated before proceeding. If a score is not changed by
  #some local move, recompute is FALSE and it is not computed on the next step
  moves$recompute <- TRUE
  moves$score <- NA
  
  return(moves)
}

#' Initialize a graph with whitelisted arcs
#'
#' @inheritParams initializeMoves
#'
#' @return bn object with specified nodes and whitelisted arcs, otherwise empty
#' 
#' @export
#'
#' @examples
initializeDag <- function(node.names, whitelist){
  dag <- empty.graph(node.names)
  if(!is.null(whitelist)){
    for(i in 1:nrow(whitelist))
    {
      dag <- set.arc(dag.current, 
                     as.character(whitelist[i,1]), 
                     as.character(whitelist[i,2]))
    }
  }
  return(dag)
}

#' Initialize cluster for DAG fitting
#'
#' Starts parallel cluster and loads bnlearn library, functions in Scoring and
#' Fit DAG files and data
#'
#' @param dat dataframe to load to clusters
#'
#' @return parallel cluster with all necessary information preloaded
#' @export
#'
#' @examples
initializeCluster <- function(dat){
  cl <- makeCluster(detectCores() - 1)
  
  invisible(clusterEvalQ(cl, library(bnlearn)))
  invisible(clusterEvalQ(cl, source("Scoring.R")))
  invisible(clusterEvalQ(cl, source("Fit DAG.R")))
  invisible(clusterExport(cl, "dat", envir = environment()))
  return(cl)
}

#' Executes move and updates scores
#'
#' Applies a given move to a DAG and updates the node and move scores
#' (exploiting score decomposability)
#'
#' @param move character vector containing to and from node (in that order)
#' @param dag DAG to apply move to
#' @param dat dataset to use
#' @param moves dataframe of moves with scores up to date for pre-move DAG
#' @param node.scores vector of scores per node
#' @param penalty string specifying the growth rate of the complexity penalty.
#'   Either "bic", "aic" or a number between 0 and 1
#' @param cl cluster to use if computations should be done in parallel
#' @param recomputeMoves logical indicating whether scores should be recomputed
#'   If false, scores are marked for recomputing but computation is not done
#'
#' @return Returns the new dag and the new moves dataframe. If scores are
#'   recomputed, also returns the score for each node and the number of queries
#' @export
#'
#' @examples
updateHC <- function(move, dag, dat, moves, node.scores, penalty, 
                     cl = NULL, recomputeMoves = TRUE){
  node.from <- move[1]
  node.to <- move[2]
  
  #Update DAG
  moveType <- findMoveType(node.from, node.to, dag)
  dag.postMove <- makeMove(node.from, node.to, dag, moveType)
  
  #Mark nodes whose parental set has changed for updating
  nodesToRescore <- ifelse(moveType == "reverse",
                           c(node.to, node.from), c(node.to))
  
  #Mark moves operating on arcs towards nodes whose parental sets have 
  #changed for updating
  moves$recompute[moves$to == node.to] <- TRUE
  if(moveType == "reverse")
  {moves$recompute[moves$to == node.from] <- TRUE}
  
  #Mark this move for updating because its meaning has changed (e.g. if
  #move was adding an arc, it is now an arc reversal)
  moves$recompute[moves$from == node.to & moves$to == node.from] <- TRUE
  
  #Some operations that lead to a Markov equivalent graph might no longer do
  #so after the current operation was performed. (e.g. an edge reversal might
  #now lead to a new v-structure, while it previously did not). This statement 
  #ensures that moves that used to lead to Markov equivalent graphs are 
  #checked. The converse effect (move now leading to equivalent graph) is less
  #problematic, as it will just lead to a graph with the same score. It will
  #not cycle as the reverse move will be recomputed after that iteration and
  #its score difference will be valued at 0
  moves$recompute[(moves$from == node.to | moves$from == node.from) &
                    abs(moves$score) < 1E-10] <- TRUE
  
  #Update node scores and move score differences if required
  if(recomputeMoves){
    result <- recompute(dag = dag.postMove,
                        dat = dat,
                        penalty = penalty,
                        moves = moves,
                        node.scores = node.scores,
                        nodesToRescore = nodesToRescore,
                        cl = cl)
    return(list(dag = dag.postMove, moves = result$moves, 
                nodeScores = result$nodeScores, queries = result$queries))
  }
  
  return(list(dag = dag.postMove, moves = moves))
  
}

#' Updates local scores and score differences marked for updating
#'
#' Recomputes score per node for a specified subset of nodes and score
#' difference due to moves for a marked subset of moves. If no moves dataframe
#' and node scores are provided, they are initialized and computed from scratch
#'
#' @param nodesToRescore character vector specifying which nodes should have
#'   their scores updated. If not provided, scores are updated for all nodes
#' @inheritParams updateHC
#' @return Updated moves, node scores and number of queries used in computing
#'   them
#' @export
#'
#' @examples
recompute <- function(dag, dat, penalty, moves = NULL, node.scores = NULL, 
                      nodesToRescore = nodes(dag), cl, blacklist, whitelist){
  #If moves and node scores have not been provided, initialize them
  if(is.null(moves)){
    moves <- initializeMoves(nodes(dag), blacklist, whitelist)
  }
  if(is.null(node.scores)){
    node.scores <- vector("numeric", length(nodes(dag)))
    names(node.scores) <- nodes(dag)
  }
  
  #Convert moves dataframe to matrix to ensure no whitespaces are added
  #when apply does the conversion internally
  moveMat <- sapply(moves, format, trim = TRUE)
  
  #If no cluster is provided, recompute scores sequentially
  if(is.null(cl)){
    #Compute score per node
    newScores <- sapply(nodesToRescore, computeScore.node, dag = dag, 
                        dat = dat, penalty = penalty)
    node.scores[nodesToRescore] <- newScores
    
    #Update score difference for each move marked for recomputation.
    moves$score[moves$recompute] <- apply(moveMat[moves$recompute,], 1, 
                                          evaluateScore, 
                                          dat = dat, 
                                          node.scores = node.scores, 
                                          dag = dag, 
                                          penalty = penalty)
  }
  #If a cluster is provided, do same as above but in parallel
  else{
    newScores <- parSapply(cl, nodesToRescore, computeScore.node, dag = dag, 
                           dat = dat, penalty = penalty)
    node.scores[nodesToRescore] <- newScores
    
    moves$score[moves$recompute] <- parApply(cl, moveMat[moves$recompute,], 1,
                                             evaluateScore, 
                                             dat = dat, 
                                             node.scores = node.scores, 
                                             dag = dag, 
                                             penalty = penalty)
  }
  
  
  no.queries <- length(nodesToRescore) + computeNoQueries(moves, dag)
  
  #Mark all moves as updated
  moves$recompute <- FALSE
  
  return(list(moves = moves, nodeScores = node.scores, queries = no.queries))
}

#Computes number of queries involved in rescoring the moves dataframe for a 
#given DAG. Each arc addition or deletion can be scored in one query, each 
#reversal can be scored in two
computeNoQueries <- function(moves, dag){
  return(sum(apply(moves[moves$recompute,], 1, function(move){
    ifelse(findMoveType(move[1], move[2], dag) == "reverse", 2, 1)})))
}

#' Computes score change due to a move
#'
#' Compares the pre- and post-move score for nodes of which the parental set
#' changes due to the move. For arc addition or deletion this is the to node,
#' for arc reversal this is the to and the from node.
#'
#' @param move character vector containing to and from node for move
#' @param node.scores pre-move node scores
#' @inheritParams evaluateHC
#'
#' @return score difference due to move if move is valid and data is available
#'   for each parental configuration. If new graph is cyclic, returns -Inf. 
#'   If there is no data available for some parental configuration, returns NaN
#'   and prints a warning
#' @export
#' 
evaluateScore <- function(move, dat, node.scores, dag, penalty){
  #Obtain move type and DAG after applying move
  node.from <- move[1]
  if(!node.from %in% nodes(dag)) {node.from <- str_trim(node.from)}
  node.to <- move[2]
  if(!node.to %in% nodes(dag)) {node.to <- str_trim(node.to)}
  
  moveType <- findMoveType(node.from, node.to, dag)
  dag.postMove <- makeMove(node.from, node.to, dag, moveType)
  
  #If move is illegal or leads to cycle in graph, score with negative infinity
  if(class(dag.postMove) == "try-error"){return(-Inf)}
  
  #Compute difference in score from changing parent set of to-node
  newScore.to <- computeScore.node(node.to, dag = dag.postMove, dat = dat, 
                                   penalty = penalty)
  scoreDiff <- newScore.to - node.scores[[node.to]]
  
  #If the move is a reversal, also take into account change at from-node
  if(moveType == "reverse"){
    newScore.from <- computeScore.node(node.from, dag = dag.postMove, dat = dat, 
                                       penalty = penalty)
    scoreDiff <- scoreDiff + newScore.from - node.scores[node.from]
  }
  
  return(scoreDiff)
}

#' Applies move to a DAG and returns new graph
#'
#' Whether the move is an addition, deletion or reversal depends on whether the
#' arc or its reverse is already present in the DAG
#'
#' @param node.from string specifying first node in move
#' @param node.to string specifying second node in move
#' @param moveType optional string specifying whether the move is an arc
#'   addition ("add"), deletion ("delete") or reversal ("reverse")
#'
#' @return If move gives a DAG, it is returned, otherwise a "try-error" is
#'   returned
#' @export
#'
#' @examples
makeMove <- function(node.from, node.to, dag, 
                     moveType = findMoveType(node.from, node.to, dag)){
  if(moveType == "add")
  {try(set.arc(dag, node.from, node.to), silent = TRUE) %>% return(.)}
  else if(moveType == "delete")
  {drop.arc(dag, node.from, node.to) %>% return(.)}
  else if(moveType == "reverse")
  {try(reverse.arc(dag, node.from, node.to), silent = TRUE) %>% return(.)}
}


#' Determines whether a move is an arc addition, deletion or reversal for a
#' given DAG
#'
#' @inheritParams makeMove
#'
#' @return If no arc is present from node.from to node.to, returns "add". If an
#'   arc is present, returns "delete". If an arc from node.to to node.from is
#'   present, returns "reverse"
#' @export
#'
#' @examples
findMoveType <- function(node.from, node.to, dag){
  if(node.from %in% parents(dag, node.to)){return("delete")}
  else if(node.to %in% parents(dag, node.from)){return("reverse")}
  else {return("add")}
}

#' Determines whether a given move leads to a valid DAG
#'
#' Checks whether a move does not lead to directed cycles. If a tabu list is
#' provided, it checks whether the new graph is not Markov equivalent to any
#' graph in the tabu list. If a dataset over nodes is provided, it checks
#' whether data is available for each parental configuration in the new graph.
#' 
#' The last is not relevant during hill climbing or tabu search, as such moves
#' will have score -Inf, but it is relevant during random restarts to ensure 
#' that we do not start with a graph with NaN score because data is unavailable
#' for some parental configuration for a node.
#'
#' @param tabuList optional tabu list of DAGs
#' @param dat optional dataframe to be used when checking whether data is
#'   available for each parental configuration
#' @inheritParams evaluateScore
#'
#' @return TRUE if new graph meets conditions, FALSE otherwise
#' 
#' @export
#'
#' @examples
checkIfAllowed <- function(move, dag, tabuList = list(), dat = NULL){
  node.from <- move[1]
  node.to <- move[2]
  moveType <- findMoveType(node.from, node.to, dag)
  
  proposedDag <- makeMove(node.from, node.to, dag, moveType)
  
  #Check that the proposed graph is acyclic
  if(class(proposedDag) == "try-error"){return(FALSE)}
  
  #Check that data is available for all of the configurations of the proposed
  #parental set. Otherwise, the likelikhood is undefined  
  if(!is.null(dat)){
    proposedScore <- computeNAL.discrete(node.to, proposedDag, dat)
    if(length(proposedScore) == 1 && is.nan(proposedScore)){return(FALSE)}
  }
  
  #Check that proposed graph is not Markov equivalent to any graph in the
  #tabu list
  return(!graphInList(proposedDag, tabuList))
}


#' Finds allowed move with largest associated score change
#'
#' @param moves dataframe of moves with updated scores
#' @inheritParams checkIfAllowed 
#'
#' @return allowed move with largest associated score change
#' 
#' @export
#'
#' @examples
findBestAllowedMove <- function(moves, dag, tabuList = list()){
  if(any(moves$recompute))
  {warning("Not all scores are up-to-date: found move potentially suboptimal")}
  
  #Order score differences from largest to smallest
  scoreOrder <- order(moves$score, decreasing = TRUE)
  
  #Iterate from largest to smallest score and return first allowed move
  for(moveIndex in scoreOrder)
  {
    currentMove <- as.character(moves[moveIndex, ])
    
    if(checkIfAllowed(currentMove, dag, tabuList))
    {return(currentMove)}
  }
  
  #This shouldn't occur for sane tabu list lengths
  stop("No allowed moves found") 
}

#' Checks whether a DAG is Markov equivalent to any graph in a given list
#'
#' @param dag DAG to compare to graphs in list
#' @param graphs list of graphs
#'
#' @return TRUE if given DAG is Markov equivalent to any graph in list, FALSE
#'   otherwise
#' @export
#'
#' @examples
graphInList <- function(dag, graphs){
  #Check whether two DAGs are in the same Markov equivalence class, allowing
  #for either to be null (in which case the comparison returns false)
  compareDAGs <- function(dag1 , dag2){
    if(is.null(dag1)||is.null(dag2)){return(FALSE)}
    else {return(all.equal(cpdag(dag1), cpdag(dag2)) == TRUE)}
  }
  
  #If list is empty, return false by default
  if(length(graphs) < 1){return(FALSE)}
  else{
    comparison <- sapply(graphs, compareDAGs, dag2 = dag)
    return(any(comparison))
  }
}

#' Performs tabu search from some initial graph
#'
#' Given some initial graph, makes at most tabuSteps moves to graphs that are
#' not in a tabu list of length at most tabuLength. Each move is the best
#' allowed move. If any of the obtained graphs has strictly higher score than
#' the initial graph, the algorithm terminates and the new graph is returned
#'
#' If no better graph is found, NULL is returned in place of a graph. Any
#' function calling this routine should cope gracefully with this scenario
#'
#' @param tabuSteps maximum number of tabu search steps to perform
#' @param tabuLength maximum number of previously seen graphs to store
#' @param dag.start initial DAG
#' @inheritParams evaluateHC
#'
#' @return If a better DAG than the initial one is found, returns a list
#'   containing the DAG, an updated moves dataframe and node scores vector and
#'   the number of queries performed during the search procedure. If no better
#'   graph is found, returns a list containing a NULL graph and the number of
#'   queries
#'
#' @export
#'
#' @examples
tabuSearch <- function(tabuSteps, tabuLength, dag.start, dat, moves, 
                       node.scores, penalty, cl){
  tabuList <- vector(mode = "list", length = tabuLength)
  no.queries.local <- 0
  dag <- dag.start
  
  #Record score of starting graph
  bestScore <- sum(node.scores)
  
  for(tabuStep in 1:tabuSteps){
    #Add previous graph to tabu list
    tabuList[[(tabuStep - 1)%%tabuLength + 1]] <- dag
    
    #Find best scoring neighboring DAG that is not in tabu list
    bestMove <- findBestAllowedMove(moves, dag, tabuList)
    
    #Execute best move
    postMove <- updateHC(bestMove,
                         dag = dag,
                         dat = dat, 
                         moves = moves, 
                         node.scores = node.scores,
                         penalty = penalty, 
                         cl = cl)
    dag <- postMove$dag
    moves <- postMove$moves
    node.scores <- postMove$nodeScores
    no.queries.local <- no.queries.local + postMove$queries
    
    #If score of new graph exceeds score of starting graph,
    #return new graph 
    newScore <- sum(node.scores)
    if(newScore > bestScore)
    {return(list(dag = dag, moves = moves, nodeScores = node.scores,
                 queries = no.queries.local))}
  }
  #If no improvement has been found, signal this by returning
  #null and the number of queries
  return(list(dag = NULL, queries = no.queries.local))
}

#' Performs a series of random moves on a given DAG
#'
#' Given a DAG, performs a given number of random moves. Each random move is
#' checked to give a valid DAG with data available for each parental
#' configuration. After performing random moves, recomputes those scores that
#' need updating.
#'
#' As the initial DAG is often sparse and the algorithm randomly generates pairs
#' of nodes to determine the moves, it tends to favor arc additions. I am not
#' sure whether I am happy with this, I might change this later and see if it
#' gives a better performance
#'
#' @param randomMoves number of random moves to do
#' @inheritParams updateHC
#'
#' @return list containing DAG after random moves, updated moves dataframe and
#' node score vector and number of queries used for score recomputation
#' 
#' @export
#'
#' @examples
perturbGraph <- function(randomMoves, dag, dat, moves, penalty, cl){
  node.names <- nodes(dag)
  no.moves <- 0
  dag.current <- dag
  
  while(no.moves < randomMoves){
    #Generate random move
    move <- sample(node.names, 2, replace = FALSE)
    
    #Check that move gives a valid graph
    if(checkIfAllowed(move, dag.current, dat = dat)){
      no.moves <- no.moves + 1
      
      #Adapt graph and record which moves will need to be recomputed,
      #but do not recompute at each iteration (do this once at the end)
      postMove <- updateHC(move,
                           dag = dag.current, 
                           dat = dat, 
                           moves = moves, 
                           node.scores = NULL,
                           penalty = penalty, 
                           recomputeMoves = FALSE)
      dag.current <- postMove$dag
      moves <- postMove$moves
    }
  }
  
  recomputed <- recompute(dag.current, dat, penalty, moves, cl = cl )
  moves <- recomputed$moves
  node.scores <- recomputed$nodeScores
  no.queries <- recomputed$queries
  
  return(list(dag = dag.current, moves = moves, nodeScores = node.scores,
              queries = no.queries))
}



#' Learns BN structure from data using tabu search with random restarts
#'
#' Starts from the empty graph. Uses a hill climber with arc addition, deletion
#' and reversal operations to find locally optimal graph structures. Score
#' function is the likelihood minus a term penalizing graph complexity. Uses
#' tabu search and random restart meta heuristics to decrease the probability of
#' ending up in a local optimum that is far from the global optimum in score.
#' Exploits score decomposability to efficiently evaluate the score difference
#' of each neighbouring graph. Can do most costly computations in parallel
#' (works well for large n, as computing scores becomes costly in this case).
#'
#' Records number of queries involved in computation. Each query involves
#' computing the score for a single node.
#'
#' @param dat dataset to learn BN from
#' @param penalty determines penalization term on score function. Must be "bic",
#'   "aic" or a number larger than 0
#' @param parallel if TRUE, computations are done in parallel where possible,
#'   using all but one of the available threads on the machine
#' @param tabuSteps maximum number of steps to do in tabu search
#' @param tabuLength maximum length of tabu list
#' @param restarts number of random restarts to do
#' @param randomMoves number of moves to do on DAG when doing random restart
#' @inheritParams initializeMoves
#'
#' @return list containing best DAG found and number of queries involved in the
#'   computation
#' @export
#'
#' @examples
#' fit.dag(learning.test, penalty = "bic")
fit.dag <- function(dat, penalty, parallel = TRUE, cl = NULL,
                    tabuSteps = 10, tabuLength = tabuSteps,
                    restarts = 0, randomMoves = ncol(dat),
                    blacklist = NULL, whitelist = NULL,
                    dag.start = NULL){
  
  #If no graph is provided, initialize to empty graph + whitelisted edges
  if(is.null(dag.start)){
    dag.current <- initializeDag(names(dat), whitelist)
    no.queries <- 0
  }
  else{
    dag.current <- dag.start
    no.queries <- dag.start$learning$ntests
  }
  
  #Initialize cluster and load relevant data, or set cluster to null
  if(parallel && is.null(cl)){
    cl <- initializeCluster(dat)
    on.exit(stopCluster(cl))
  }
  
  #Initialize dataframe of moves and compute current score per node and change
  #in score from each move
  initialization <- recompute(dag.current, dat, penalty, cl = cl, 
                              blacklist = blacklist, whitelist = whitelist)
  moves <- initialization$moves
  node.scores <- initialization$nodeScores
  no.queries <- no.queries + initialization$queries
  
  restart <- 0
  
  #Track which DAG was the best over all restarts and what its score was
  dag.best <- dag.current
  score.best <- sum(node.scores)
  
  while(restart <= restarts){
    converged <- FALSE
    
    while(!converged){
      
      bestMove <- findBestAllowedMove(moves, dag.current)
      
      if(as.numeric(bestMove[4]) > 0){
        #Execute best move. Update DAG, scores for nodes, score differences 
        #for moves and number of contingency tables computed
        postMove <- updateHC(bestMove,
                             dag = dag.current, 
                             dat = dat, 
                             moves = moves, 
                             node.scores = node.scores,
                             penalty = penalty, 
                             cl = cl)
        dag.current <- postMove$dag
        moves <- postMove$moves
        node.scores <- postMove$nodeScores
        no.queries <- no.queries + postMove$queries
      }
      else if (tabuSteps > 0){
        tabuResult <- tabuSearch(tabuSteps = tabuSteps,
                                 tabuLength = tabuLength,
                                 dag.start = dag.current, 
                                 dat = dat, 
                                 moves = moves, 
                                 node.scores = node.scores,
                                 penalty = penalty,
                                 cl = cl)
        
        #If tabu search returned null, it did not find a DAG with higher
        #score than the previous best. In this case the search terminates 
        #(or does a random restart) from the previous best graph. 
        #Otherwise, we continue searching from the new best graph.
        if(is.null(tabuResult$dag)){converged <- TRUE}
        else {
          dag.current <- tabuResult$dag
          moves <- tabuResult$moves
          node.scores <- tabuResult$nodeScores
        }
        no.queries <- no.queries + tabuResult$queries
      }
      #If no score-increasing move exists and tabu search is not active,
      #terminate current phase of search
      else{converged <- TRUE}
    }
    
    score.current <- computeScore(dag.current, dat, penalty, cl = cl)
    if(score.current > score.best){
      dag.best <- dag.current
      score.best <- score.current
      
      #If a better graph is found, reset number of restarts
      restart <- 0
    }
    
    #Only go through perturbation procedure if max number of restarts
    #has not been reached yet.
    if(restart < restarts){
      #Perturb graph
      perturbed <- perturbGraph(randomMoves,
                                dag = dag.best,
                                dat = dat,
                                moves = moves,
                                penalty = penalty,
                                cl = cl)
      
      #Retrieve perturbed graph and updated moves and node scores
      dag.current <- perturbed$dag
      moves <- perturbed$moves
      node.scores <- perturbed$nodeScores
      no.queries <- no.queries + perturbed$queries
    }
    
    restart <- restart + 1
  }
  
  dag.best$learning$ntests <- no.queries
  return(dag = dag.best)
}


#' Fits DAG with known node order and restricted size of parental sets
#'
#' @param node.names character vector of ordered node names
#' @param max.parents maximum number of parents per node
#' @inheritParams fit.dag
#'
#' @return Maximum scoring DAG for given restrictions on node order and parental
#'   sets
#'
#' @export
#'
#' @examples
fit.dag.ordered <- function(node.names, max.parents, dat, penalty, 
                            parallel = TRUE){
  
  #Initialize DAG and cluster
  dag.fitted <- empty.graph(node.names)
  no.nodes <- length(node.names)
  
  if(parallel){
    cl <- initializeCluster(dat)
    on.exit(stopCluster(cl))
  }
  
  #Computes score for given parent configuration of a node. Returns -Inf if
  #there is no data for some parental configuration
  wrapperScore <- function(parents.node, node, dat, penalty, no.nodes){
    
    dag <- empty.graph(c(node, parents.node))
    if(length(parents.node) > 0) {parents(dag, node) <- parents.node}
    
    score.node <- computeScore.node(node, dag, dat, penalty, 
                                    no.nodes = no.nodes)
    
    #If any parent configuration has no occurences, we deem it to be 
    #too complex to learn from our dataset, and give it -Inf score
    if(is.na(score.node)){score.node <- -Inf}
    
    return(score.node)
  }
  
  #Determine optimal parental sets for each node separately (possible due to
  #known node ordering, which guarantees acyclicity when only nodes earlier in
  #the ordering are considered as candidate parents)
  for(i in 2:no.nodes){
    #Generate list of potential parental sets
    node <- node.names[i]
    possible.parents <- lapply(0:min(max.parents, i-1), 
                               combn, x = node.names[1:(i-1)],
                               simplify = FALSE) %>% 
      unlist(., recursive = FALSE)
    
    #Compute score for each possible parental configuration.
    #Only worth incurring paralellization overhead if set of possible
    #parents is large enough 
    if(parallel && length(possible.parents) > 100){
      scores <- parSapply(cl, possible.parents, 
                          FUN = wrapperScore,
                          node = node,
                          dat = dat,
                          penalty = penalty,
                          no.nodes = length(node.names))
    }
    else{
      scores <- sapply(possible.parents,
                       FUN = wrapperScore,
                       node = node,
                       dat = dat,
                       penalty = penalty,
                       no.nodes = length(node.names))
    }
    
    bestParents <- possible.parents[[which.max(scores)]]
    
    if(length(bestParents) > 0)
    {parents(dag.fitted, node) <- bestParents}
  }
  
  #Compute and store number of parental configurations that needed to be
  #recorded
  no.queries <- sum(sapply(0:max.parents, function(x){choose(1:no.nodes, x)}))
  dag.fitted$learning$ntests <- no.queries
  return(dag.fitted)
}


#' Impute data in parallel using a fitted Bayesian Network
#'
#' Distributes a dataframe with missing values over all but one of the clusters
#' on the current machine and uses likelihood weighted sampling from the
#' mutilated network to impute missing values by their posterior mean
#' (continuous) or mode (discrete).
#'
#' @param dag.fit Fitted Bayesian Network to use for data imputation
#' @param dat Dataframe with missing values to impute
#' @param particles Number of samples to use in simulation
#'
#' @return dat with missing values imputed
#' @export
#'
#' @examples
parImpute <- function(dag.fit, dat, particles){
  #Initialize parallel cluster
  no.cores <- detectCores() - 1
  cl <- makeCluster(no.cores)
  on.exit(stopCluster(cl))
  invisible(clusterEvalQ(cl, library(bnlearn)))
  
  #Split dataframe into equal sized chunks
  rownames(dat) <- 1:nrow(dat)
  no.obs <- ceiling(nrow(dat)/no.cores)
  dat.list <- split(dat, (as.numeric(rownames(dat))-1)%/%no.obs)
  
  #Impute in parallel and combine and return result
  impute.list <- parLapply(cl, dat.list, impute, object = dag.fit, 
                           method = "bayes-lw", n = particles)
  imputed <- rbindlist(impute.list)
  return(imputed)
}


#' Performs parametric EM for a given DAG
#'
#' Initializes by learning parameters using locally complete data, then imputes
#' missing values (E-step) and refits parameters (M-step) a given number of
#' times using likelihood weighted sampling. Allows setting a fixed maximum
#' number of observations per parameter for computational reasons. These
#' observations are randomly sampled from the full dataset.
#'
#' @param dag DAG structure to learn parameters for
#' @param dat Dataframe with missing values to learn parameters from
#' @param max.iter Maximum number of iterations (default is 1)
#' @param particles Number of samples to use in simulation
#' @param obs.perParam Number of observations per parameter to learn. Defaults
#'   to full dataset
#' @param parallel If TRUE, do data imputation in parallel
#'
#' @return DAG with fitted parameters and imputed data after last E-step
#'
#' @export
#'
#' @examples
em.parametric <- function(dag, dat, max.iter = 1, particles = 500, 
                          obs.perParam = NULL, parallel = FALSE){
  #Initialize by fitting parameters using locally complete data
  if(class(dag)[1] != "bn.fit"){
    if(class(dag)[1] != "bn"){stop("DAG should be either of type bn or type bn.fit")}
    dag.fit <- bn.fit(dag, dat, replace.unidentifiable = TRUE)
  }
  else{
    dag.fit <- dag
  }
  dat.imputed <- NULL
  converged <- FALSE
  iteration <- 1
  
  #Take subsample of observations if desired
  if(!is.null(obs.perParam)){
    no.obs <- min(obs.perParam*nparams(dag.fit), nrow(dat))
    dat <- sample_n(dat, size = no.obs)
  }
  
  #Iterate until maximum number of iterations is reached
  while(!converged && iteration <= max.iter){
    iteration <- iteration + 1
    
    if(parallel){
      dat.imputed <- parImpute(dag.fit, dat, particles = particles)
    } else {
      dat.imputed <- impute(dag.fit, dat, method = "bayes-lw", n = particles)
    }
    
    dag.fit <- bn.fit(dag, dat.imputed)
    #TODO: add criterion for convergence
  }
  
  return(list(dag = dag.fit, dat = dat.imputed))
}

#' Uses Structural EM to learn DAG structure from incomplete data
#'
#' Iteratively imputes data (E-step) and learns DAG structure from fitted data
#' (M-step). Continuous until a maximum number of iterations is reached or the
#' SHD between the new and the previous graph is zero.
#'
#' @param dat Incomplete dataframe to learn DAG from
#' @param max.iter Maximum number of iterations
#' @param penalty Penalization factor to use for scoring, default is "bic"
#' @inheritParams fit.dag
#'
#' @return Fitted DAG
#' @export
#'
#' @examples
em.structural <- function(dat, tabuSteps = 10, penalty = "bic",
                          tabuLength = tabuSteps, blacklist = NULL, 
                          whitelist = NULL, max.iter = Inf, particles = 500,
                          parallel = FALSE){
  #Initialize to graph with only whitelisted edges and 
  dag.current <- initializeDag(names(dat), whitelist)
  dat.imputed <- dat
  fitted.current <- bn.fit(dag.current, dat.imputed)
  
  converged <- FALSE
  iteration <- 1
  
  while(!converged && iteration <= max.iter){
    iteration <- iteration + 1
    
    #E-step: fit parameters of current network using previously imputed data
    #and impute data based on new network
    if(parallel){
      dat.imputed <- parImpute(fitted.current, dat, particles = particles)
    } else {
      dat.imputed <- impute(fitted.current, dat, method = "bayes-lw", n = particles)
    }
    
    #M-step: optimize using completed data
    dag.new <- fit.dag(dat.imputed, penalty = penalty,
                       parallel = parallel,
                       tabuSteps = tabuSteps,
                       tabuLength = tabuLength,
                       blacklist = blacklist,
                       whitelist = whitelist,
                       dag.start = dag.current)
    
    #Use Bayesian fitting if all data is discrete to avoid unidentiability.
    #If continuous data is also present, use MLE
    
    if(all(sapply(dat, is.factor))){
      fitted.new <- bn.fit(dag.new, dat.imputed, method = "bayes", iss = 1)
    }
    else{
      fitted.new <- bn.fit(dag.new, dat.imputed, replace.unidentifiable = TRUE)
    }
    
    if(shd(dag.current, dag.new) == 0){converged <- TRUE}
    else{dag.current <- dag.new; fitted.current <- fitted.new}
  }
  
  return(dag.new)
}


#' Learns DAGs using the non-parametric bootstrap
#'
#' Repeatedly sample with replacement from original dataset and learn DAGs from
#' these samples. 
#'
#' @param samples Number of DAGs to learn
#' @param dat Data to obtain bootstrap samples from
#' @param sampleSize Number of observations to sample. Default is nrow(dat)
#' @param parallel Bootstrap in parallel if TRUE (default)
#' @inheritParams fit.dag
#'
#' @return List of bn objects learned from bootstrap samples
#' @export
#'
#' @examples
bootstrap.dag <- function(samples, dat, penalty, 
                          tabuSteps = 10, tabuLength = tabuSteps, 
                          blacklist = NULL, whitelist = NULL,  
                          sampleSize = nrow(dat), parallel = FALSE){
  require(dplyr)
  
  #Obtains bootstrap sample and learns and returns DAG
  bootstrap.single <- function(dummy, dat, penalty, tabuSteps, tabuLength,
                               blacklist, whitelist, sampleSize){
    dat.sample <- sample_n(dat, sampleSize, replace = TRUE)
    dag.fitted <- fit.dag(dat.sample, penalty = penalty,
                          parallel = FALSE,
                          tabuSteps = tabuSteps,
                          tabuLength = tabuLength,
                          blacklist = blacklist,
                          whitelist = whitelist)
    return(dag.fitted)
  }
  
  #Run in parallel if requested
  if(parallel){
    
    cl <- initializeCluster(dat)
    on.exit(stopCluster(cl))
    dag.samples <- parLapplyLB(cl, 1:samples, bootstrap.single, dat = dat, 
                               penalty = penalty,
                               tabuSteps = tabuSteps,
                               tabuLength = tabuLength,
                               blacklist = blacklist,
                               whitelist = whitelist,
                               sampleSize = sampleSize) 
  } else {
    dag.samples <- lapply(1:samples, bootstrap.single, dat = dat, 
                          penalty = penalty,
                          tabuSteps = tabuSteps,
                          tabuLength = tabuLength,
                          blacklist = blacklist,
                          whitelist = whitelist,
                          sampleSize = sampleSize)
  }
  
  return(dag.samples)
}

