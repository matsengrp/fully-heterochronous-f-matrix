source("beta2.R")

beta_chain_integral <- function(a, b, alpha, beta) {
  # a, b: numeric vectors of length n, entries should be 0 or 1
  # alpha, beta: positive shape parameters of the Beta(alpha, beta)
  
  ## -------- basic checks --------
  if (length(a) != length(b)) {
    stop("a and b must have the same length.")
  }
  n <- length(a)
  if (n == 0L) {
    stop("Length of a and b must be at least 1.")
  }
  if (!all(a %in% c(0, 1))) {
    stop("All entries of a must be 0 or 1.")
  }
  if (!all(b %in% c(0, 1))) {
    stop("All entries of b must be 0 or 1.")
  }
  
  a <- as.integer(a)
  b <- as.integer(b)
  
  # Precompute tail sums of a: A_k = sum_{r=k}^n a_r
  A <- rev(cumsum(rev(a)))             # length n
  B0 <- beta(alpha, beta)
  
  # Indices where epsilon can vary (b_r = 1)
  idx_free <- which(b == 1L)
  m <- length(idx_free)
  
  # Maximum possible E_k = sum_{r=k}^n (a_r + eps_r) is bounded by sum(a + b)
  maxE <- sum(a + b)
  
  # Precompute all Beta(alpha + t, beta) / Beta(alpha, beta) for t = 0, ..., maxE
  Bt <- beta(alpha + 0:maxE, beta) / B0   # length maxE + 1, index by t+1
  
  ## -------- special case: no free eps (all b_r = 0) --------
  if (m == 0L) {
    # S_k(eps) = 0 for all k, so E_k = A_k
    # result = prod_k Beta(alpha + A_k, beta) / Beta(alpha, beta)
    return(prod(Bt[A + 1L]))
  }
  
  ## -------- generate all epsilon patterns on free coordinates --------
  # eps_free: 2^m x m matrix, each row an epsilon pattern on idx_free
  eps_free <- as.matrix(expand.grid(rep(list(0L:1L), m)))
  storage.mode(eps_free) <- "integer"
  n_eps <- nrow(eps_free)
  
  ## -------- compute S_k(eps) for all epsilon at once --------
  # S_k(eps) = sum_{r >= k} eps_r
  # Only positions in idx_free can be 1; define T_{k,j} = 1 if idx_free[j] >= k
  k_idx <- seq_len(n)
  T <- outer(k_idx, idx_free, FUN = function(k, r) as.integer(r >= k))  # n x m
  
  # S_all: n_eps x n matrix, S_all[i, k] = S_k(eps^(i))
  # (matrix multiplication uses BLAS and is much faster than R loops)
  S_all <- eps_free %*% t(T)            # (n_eps x m) %*% (m x n) = (n_eps x n)
  storage.mode(S_all) <- "integer"
  
  ## -------- build all E_k(eps) and corresponding Beta products --------
  # Ek_all[i, k] = A_k + S_k(eps^(i))
  Ek_all <- matrix(rep(A, times = n_eps), nrow = n_eps, byrow = TRUE) + S_all
  
  # Lookup Beta ratios: factor_mat[i, k] = Beta(alpha + Ek_all[i,k], beta) / B0
  factor_mat <- matrix(Bt[Ek_all + 1L], nrow = n_eps, ncol = n)
  
  # Product over k for each epsilon pattern:
  # term_prod[i] = prod_k factor_mat[i, k]
  # Use logs for numerical stability and speed
  term_prod <- exp(rowSums(log(factor_mat)))
  
  ## -------- sign for each epsilon pattern --------
  # sign(eps) = (-1)^{sum_r eps_r}
  parity <- rowSums(eps_free) %% 2L
  sign_vec <- ifelse(parity == 0L, 1, -1)
  
  ## -------- final sum over all epsilon --------
  result <- sum(sign_vec * term_prod)
  return(result)
}

beta2_likelihood <- function(mat,alpha,beta){
  
  total_prob = 1
  m = dim(mat)[1]
  n_tips = (m+2)/2
  
  for (i in 2:(m - 1)) {
    F_ii <- mat[i , i ]
    
    previous = mat[i-1,i-1]
    diagonals <- valid_heterochronous_diagonals(i, n_tips, previous)
    f_l<- diagonals[1]
    f_u <- diagonals[2]
    prob_ii = ifelse(f_l == f_u , 1, 
                     beta(alpha + (f_u - F_ii)/2, beta+(F_ii-f_l)/2 ) / beta(alpha,beta))
    total_prob = total_prob * prob_ii
  }
  
  # Fill entries off the diagonal and subdiagonal: we do one column at a time, to collect all 
  # dependent terms under the integration 
  for (c in 1:(m-3)){
    vec_a = c()
    vec_b = c()
    n=0
    for (r in (c+2):m){
      F_ij = mat[r,c] 

      left <- if (c == 1) 0 else mat[r, c - 1]
      left_up <- if (c == 1) 0 else mat[r - 1, c - 1]
      up <- mat[r - 1, c]
      prev_diag <- mat[r - 1, r - 1]
      
      non_diag <- valid_non_diagonal_values(r, c, left, left_up, up, prev_diag, isochronous = FALSE)
      f_l <- non_diag[1]
      f_u <- non_diag[2]
      
      if (f_l == f_u){next}
      n=n+1 
      vec_a = c(vec_a, f_u - F_ij)
      vec_b = c(vec_b,F_ij-f_l)
    }
    
    #just to see if it gives me the right thing
    #print(vec_a)
    
    if (n == 0) {next}
    prob_col = beta_chain_integral(vec_a, vec_b, alpha, beta)
    total_prob = total_prob*prob_col
  }
  
  
  # for (i in 3:m) {
  #   for (j in 1:(i - 2)) {
  #     F_ij = mat[i,j] 
  #     
  #     left <- if (j == 1) 0 else mat[i, j - 1]
  #     left_up <- if (j == 1) 0 else mat[i - 1, j - 1]
  #     up <- mat[i - 1, j]
  #     prev_diag <- mat[i - 1, i - 1]
  #     
  #     non_diag <- valid_non_diagonal_values(i, j, left, left_up, up, prev_diag, isochronous = FALSE)
  #     f_l <- non_diag[1]
  #     f_u <- non_diag[2]
  #     prob_ij = ifelse(f_l == f_u , 1, 
  #                      beta(alpha + f_u - F_ij, beta+F_ij-f_l) / beta(alpha,beta))
  #     total_prob = total_prob * prob_ij
  #     
  #   }
  # }
  
  return (total_prob)
}

test=matrix(c(2,0,0,0,0,0,1,3,0,0,0,0,0,2,4,0,0,0,0,2,3,3,0,0,0,2,2,2,2,0,0,1,1,1,1,1),nrow=6, ncol=6,byrow=TRUE)
beta2_likelihood(test,0.3,0.7)




source("~/Documents/Phd/25 Fall/farmrats_working/Rscript/enumerate_clean.R")

alpha = 1.1
beta=1.5
prob= rep(0,length(F.list))
for (i in 1:length(F.list)){
  mat = F.list[[i]]
  prob[i]=beta2_likelihood(mat,alpha,beta )
}

#out = construct_trees_het_D(20)


# beta2_min_max_mat <-function(alpha,beta,n){
#   
# }

