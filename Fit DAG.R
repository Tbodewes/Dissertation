library(bnlearn)
library(parallel)
library(stringr)
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
#' @param node.names character vectors of nodes
#'
#' @return dataframe containing each move, its score and whether the score is
#'   up-to-date
#' @export
#'
#' @examples
initializeMoves <- function(node.names){
  #Get all pairs of distinct nodes
  moves <- expand.grid(from = node.names, to = node.names, stringsAsFactors = F)
  moves <- with(moves, moves[from != to,])
  
  #Indicates that none of the scores for the moves has been computed, hence
  #all should be evaluated before proceeding. If a score is not changed by
  #some local move, recompute is FALSE and it is not computed on the next step
  moves$recompute <- TRUE
  moves$score <- NA
  
  return(moves)
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
                      nodesToRescore = nodes(dag), cl){
  #If moves and node scores have not been provided, initialize them
  if(is.null(moves)){
    moves <- initializeMoves(nodes(dag))
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
    newCounts <- lapply(nodesToRescore, computeCount.node, dag = dag, 
                        dat = dat)
    newScores <- sapply(newCounts, computeScore.node, penalty = penalty,
                        no.nodes = nnodes(dag))
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
    newCounts <- parLapply(cl, nodesToRescore, computeCount.node, dag = dag, 
                           dat = dat)
    newScores <- parSapply(cl, newCounts, computeScore.node, 
                           penalty = penalty, no.nodes = nnodes(dag))
    node.scores[nodesToRescore] <- newScores
    
    moves$score[moves$recompute] <- parApply(cl, moveMat[moves$recompute], 1,
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

#Computes number of queries involved in rescoring the moves dataframe. Each
#arc addition or deletion can be scored in one query, each reversal in two
computeNoQueries <- function(moves, dag){
  return(sum(apply(moves[moves$recompute,], 1, function(move){
    ifelse(findMoveType(move[1], move[2], dag) == "reverse", 2, 1)})))
}

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
  count.to <- computeCount.node(node.to, dag.postMove, dat = dat)
  newScore.to <- computeScore.node(count.to, penalty, no.nodes = nnodes(dag))
  scoreDiff <- newScore.to - node.scores[[node.to]]
  
  #If the move is a reversal, also take into account change at from-node
  if(moveType == "reverse"){
    count.from <- computeCount.node(node.from, dag.postMove, dat = dat)
    newScore.from <- computeScore.node(count.from, penalty, no.nodes = nnodes(dag))
    scoreDiff <- scoreDiff + newScore.from - node.scores[node.from]
  }
  
  return(scoreDiff)
}

makeMove <- function(node.from, node.to, dag, moveType){
  if(moveType == "add")
  {try(set.arc(dag, node.from, node.to), silent = TRUE) %>% return(.)}
  else if(moveType == "delete")
  {drop.arc(dag, node.from, node.to) %>% return(.)}
  else if(moveType == "reverse")
  {try(reverse.arc(dag, node.from, node.to), silent = TRUE) %>% return(.)}
}

findMoveType <- function(node.from, node.to, dag){
  if(node.from %in% parents(dag, node.to)){return("delete")}
  else if(node.to %in% parents(dag, node.from)){return("reverse")}
  else {return("add")}
}

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
    proposedCount <- computeCount.node(node.to, proposedDag, dat)
    proposedScore <- computeNAL.node(proposedCount)
    if(is.nan(proposedScore)){return(FALSE)}
  }
  
  #Check that proposed graph is not Markov equivalent to any graph in the
  #tabu list
  return(!graphInList(proposedDag, tabuList))
}

findBestAllowedMove <- function(moves, dag, tabuList = list()){
  scoreOrder <- order(moves$score, decreasing = TRUE)
  
  bestMove <- NULL
  
  for(moveIndex in scoreOrder)
  {
    currentMove <- as.character(moves[moveIndex, ])
    
    if(checkIfAllowed(currentMove, dag, tabuList))
    {bestMove <- currentMove; break }
  }
  
  return(bestMove)
}

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

perturbGraph <- function(randomMoves, dag, dat, moves, penalty, cl){
  node.names <- nodes(dag)
  no.moves <- 0
  dag.current <- dag
  while(no.moves < randomMoves){
    move <- sample(node.names, 2, replace = FALSE)
    
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
  
  return(list(dag = dag.current, moves = moves, nodeScores = node.scores))
}



fit.dag <- function(dat, penalty, parallel = TRUE, 
                    tabuSteps = 10, tabuLength = tabuSteps,
                    restarts = 0, randomMoves = ncol(dat)){
  #Initialize current graph to empty graph
  dag.current <- empty.graph(names(dat))
  
  #Initialize cluster and load relevant data, or set cluster to null
  if(parallel){
    cl <- initializeCluster(dat)
    on.exit(stopCluster(cl))
  }
  else{cl <- NULL}
  
  #Initialize dataframe of moves and compute current score per node and change
  #in score from each move
  initialization <- recompute(dag.current, dat, penalty, cl = cl)
  moves <- initialization$moves
  node.scores <- initialization$nodeScores
  no.queries <- initialization$queries
  
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

    #Perturb graph
    perturbed <- perturbGraph(randomMoves,
                              dag = dag.best,
                              dat = dat,
                              moves = moves,
                              penalty = penalty,
                              cl = cl)
    
    
    dag.current <- perturbed$dag
    moves <- perturbed$moves
    node.scores <- perturbed$nodeScores
    if(any(is.nan(node.scores))){browser()}
    
    restart <- restart + 1
  }
  
  return(list(dag = dag.best, queries = no.queries))
}
