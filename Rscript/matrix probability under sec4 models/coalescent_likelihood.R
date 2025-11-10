coal_likelihood <- function(mat){
  m=dim(mat)[1]
  n = m/2 +1 
  prob = 1
  
  # after first two leaf-leaf merge
  n_leaves = 2 * n -  4 
  n_ranked_nodes = 2
  
  for (i in 3:(2*n-2)){
    denom = choose(n_leaves,2) + choose(n_ranked_nodes,2)
    cur_col = (2*n-1)-i 
    prev_diag = mat [2*n-i ,2*n-i]
    cur_diag = mat [2*n-1-i ,2*n-1-i]
    
    if (prev_diag == cur_diag-1){ # if...: this is when two leaves are merged
      prob = prob * choose(n_leaves,2)/ denom
      n_leaves = n_leaves-2
      n_ranked_nodes = n_ranked_nodes+1
      print(paste0("leaf merge:",prob,"; n_leaves:", n_leaves, "; n_nodes:", n_ranked_nodes))
      
    } else{# else: this is when two ranked nodes are merged
      prob = prob/denom
      n_ranked_nodes = n_ranked_nodes-1
      print(paste0("node merge:", prob,"; n_leaves:", n_leaves, "; n_nodes:", n_ranked_nodes))
    }
    
  }
  
  return( prob )
}


test= matrix(c(2,0,0,0,1,1,0,0,0,0,2,0,0,0,1,1),nrow=4, byrow=TRUE)
coal_likelihood(test)

test2= matrix( c(2,0,0,0,0,0,1,3,0,0,0,0,0,2,4,0,0,0,0,1,3,3,0,0,0,0,2,2,2,0,0,0,1,1,1,1), 
               nrow=6, byrow=TRUE)
coal_likelihood(test2)



