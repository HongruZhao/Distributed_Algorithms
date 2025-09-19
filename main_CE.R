# -------------------------------
# get task ID from slurm
# --------------------------------
ID_index <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# Number of parallel tasks: This corresponds to #SBATCH --array=1-20
# --------------------------------
num_parallel=25

library(glmnet)

# --------------------------------------------
# Set model parameters 
# --------------------------------------------

nreps <- 100
# nreps should be exactly divisible by num_parallel.
p_vec=c(10,50)
Model_vec = c(1,2,3)
#n_vec=c(100, 300, 500, 1000, 2000)
n_vec=c(100, 300,500,1000,2000)


#params <- expand.grid("p" = rep(p_vec, each=nreps), "Model" = Model_vec, "n"=n_vec)
params <- expand.grid("p" = p_vec , "Model" = Model_vec, "n"=n_vec)

### functions we need
{
  ######### Basic functions
  
  #####Design Matrix Generation
  ## U
  design_U <- function(m) {
    U = matrix(0, nrow = m, ncol = m-1 )
    for(j in 1:m-1){
      U[1:(j+1),j]=c(rep(1,j), - j )/sqrt(j*(j+1))
    }  
    return(U)
  }
  
  design_H_k <- function(J,k) {
    q=length(J)
    for_ind=combn(q, k)
    rep_full=ncol(for_ind)
    ind_end=rep(1,rep_full)
    for(j in 1:rep_full){
      current=1
      k_vec=for_ind[,j]
      ind_end[j]=prod(J[k_vec]-1)
      if( ind_end[j]==0 ){next}
      for(i in 1:q ){
        if( i %in% k_vec ){
          current = kronecker( design_U(J[i]) , current )
        }else{
          current = kronecker( rep(1,J[i])/sqrt(J[i]) , current )
        }
      }
      if(j==1){
        H=current
      }else{H=cbind(H,current)}
    }
    return( list(H=H, ind_end=ind_end ) )
  }
  
  design_H <- function(J,d) {
    ind_end=1
    H=rep(1,prod(J))/sqrt(prod(J))
    if(d==0){
      return( list(H=t(t(H)), ind_start=1, ind_end=1 ) )
    }else{
      for(k in 1:d ){
        output=design_H_k(J, k )  
        H=cbind(H, output$H )
        ind_end=c(ind_end, output$ind_end )
      }
      ind_sta=rep(1,length(ind_end))
      for(i in 2:length(ind_end)){
        ind_end[i] = ind_end[i] +ind_end[i-1] 
        ind_sta[i] = ind_end[i-1] + 1
      }
      colnames(H)<-NULL
      return( list(H=H, ind_start=ind_sta, ind_end=ind_end ) )
    }
  }
  
  ##### Beta Generation
  beta_generate <- function(L_sum,p,ind_end,s,unif_range=c(-2,2), spectral_normalize=5){
    beta=runif(L_sum*p, min=unif_range[1], max=unif_range[2])
    dim(beta)=c(L_sum,p)
    if( s+1< length(ind_end) ){beta[(ind_end[s+1]+1):L_sum,]=0}
    beta=beta/norm(beta, type = c("2"))*spectral_normalize
    return(beta)
  }
  ##### Beta Generation: scheme 1 and 2 and 3
  ### scheme==1, Mutually independent 
  ### scheme==2, Joint independence, 1,2,..,{q-2,q-1,q} are independent
  ### scheme==3, Conditional independence, given q,  1,2,..,{q-2,q-1} are independent
  beta_generate_scheme <- function(J,d=3,L_sum,p,ind_start_H,ind_end_H,unif_range=c(-2,2),scheme=1 ){
    beta=matrix(0,nrow=L_sum, ncol=p  )
    q=length(J)
    k_ind=0
    col_ind=c(1,sort(sample(2:p,2)))
    ii=1
    if(scheme==1 ){
      for_ind=combn(q, 1)
      rep_full=ncol(for_ind)
      for(j in 1:rep_full){
        ii=ii+1
        beta[ind_start_H[ii]:ind_end_H[ii],col_ind]=matrix(runif(length(ind_start_H[ii]:ind_end_H[ii])*length(col_ind),min=unif_range[1],max=unif_range[2]),ncol=length(col_ind))
      }
    }
    if(scheme==2){
      for(k in 1:d ){
        for_ind=combn(q, k)
        rep_full=ncol(for_ind)
        for(j in 1:rep_full){
          ii=ii+1
          if( all( for_ind[,j] %in% c(q,q-1,q-2) ) | k==1 ){
            beta[ind_start_H[ii]:ind_end_H[ii],col_ind]=matrix(runif(length(ind_start_H[ii]:ind_end_H[ii])*length(col_ind),min=unif_range[1],max=unif_range[2]),ncol=length(col_ind))
          }
        }
      }
    }
    if(scheme==3){
      for(k in 1:d ){
        for_ind=combn(q, k)
        rep_full=ncol(for_ind)
        for(j in 1:rep_full){
          ii=ii+1
          if( all( for_ind[,j] %in% c(q,q-1,q-2) ) | k==1 | (k==2 & q %in% for_ind[,j] ) ){
            beta[ind_start_H[ii]:ind_end_H[ii],col_ind]=matrix(runif(length(ind_start_H[ii]:ind_end_H[ii])*length(col_ind),min=unif_range[1],max=unif_range[2]),ncol=length(col_ind))
          }
        }
      }
    }
    return(beta)
  }
  
  ##### X Generation
  cov_decay <- function(p , r) {
    # Initialize an empty n x n matrix
    matrix <- matrix(0, p, p)
    # Populate the matrix
    for (i in 1:p) {
      for (j in 1:p ) {
        matrix[i, j] <- r^abs(i - j)
      }
    }
    return(matrix)
  }
  ##### Y Generation Poisson
  poisson_Y_generate <- function(H, X, beta) {
    u_mat = H%*%beta%*%X
    Y=rpois(length(u_mat), exp(u_mat))
    dim(Y)=dim(u_mat)
    return(Y)
  }
  ##### Y Generation multinomial
  multinomial_Y_generate <- function(H, X, beta) {
    u_mat = H%*%beta%*%X
    dim(rep(1,nrow(u_mat))%*%t(colSums(exp( u_mat ))) )
    p_mat=exp( u_mat ) / (rep(1,nrow(u_mat))%*%t(colSums(exp( u_mat ))) )
    Y=rmultinom( 1, size=1, prob=p_mat[,1])
    for(j in 2:ncol(u_mat) ){
      Y=cbind(Y,rmultinom( 1, size=1, prob=p_mat[,j]))
    }
    return(Y)
  }
  ##### Log-likelihood
  ### Log-likelihood
  ##  Poisson
  ll_pois <- function(beta,H,X,Y){
    n_sa=ncol(Y)
    u_mat=H%*% (beta %*% X)
    C=sum(exp(u_mat))
    return( (-sum(Y*u_mat) + C)/n_sa)
  }
  grad_ll_pois <- function(beta,H,X,Y){
    mu_mat=exp(H%*% (beta %*% X) )
    n_sa=ncol(Y)
    dis_mat=mu_mat - Y
    term=0
    for(j in 1:n_sa){
      term=term+dis_mat[,j]%*%t(X[,j])
    }
    return(t(H)%*%term/n_sa)
  }
  ##  multinomial n_i=1 version
  ll_mult <- function(beta,H,X,Y){
    n_sa=ncol(Y)
    u_mat=H%*% (beta %*% X)
    C_vec=log(colSums(exp(u_mat)))
    C=colSums(Y)%*%C_vec
    return( (-sum(Y*u_mat) + C)/n_sa )
  }
  
  ll_mult_se <- function(beta,H,X,Y){
    n_sa=ncol(Y)
    u_mat=H%*% (beta %*% X)
    C_vec=log(colSums(exp(u_mat)))
    C=sum(colSums(Y)*C_vec)
    ll_vec=-colSums(Y*u_mat) +colSums(Y)*C_vec
    return( list(mean_ll=mean(ll_vec), se_ll=sd(ll_vec)  ) )
  }
  
  ll_mult_theta <- function(theta,X,Y){
    n_sa=ncol(Y)
    u_mat=theta %*% X 
    C_vec=log(colSums(exp(u_mat)))
    C=colSums(Y)%*%C_vec
    return( (-sum(Y*u_mat) + C)/n_sa )
  }
  
  grad_ll_mult <- function(beta,H,X,Y){
    expu_mat=exp(H%*% (beta %*% X) )
    n_sa=ncol(Y)
    #expsum_vec= colSums(exp(u_mat))
    term_mat=(rep(1,nrow(expu_mat))%*%t(colSums( expu_mat )) )
    p_mat=expu_mat/term_mat
    dis_mat=p_mat-Y
    term=0
    for(j in 1:n_sa){
      term=term+dis_mat[,j]%*%t(X[,j])
    }
    return(t(H)%*%term/n_sa)
  }
  
  ##### group lasso
  group_lasso <- function(beta, ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
    # Create a list of row and column indices
    row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
    col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
    
    # Function to extract submatrix and compute Frobenius norm
    get_norm <- function(r, c) {
      sub_matrix <- beta[r, c, drop = FALSE]
      norm(sub_matrix, type = "F")
    }
    # Apply the function over all combinations of row and column indices
    norms <- outer(row_indices, col_indices, Vectorize(get_norm))
    return( sum(norms*w) )
  }
  
  group_lasso_proximal <- function(beta_z, w, lambda, ind_start_H, ind_end_H, ind_start_X, ind_end_X) {
    # Create lists of row and column indices
    row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
    col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
    if (length(w) == 1){
      ww=rep(w,length(row_indices)*length(col_indices))
      dim(ww)=c(length(row_indices), length(col_indices))
      w=ww
    }
    # Function to process each submatrix
    process_submatrix <- function(r, c, i, j) {
      b_sub <- beta_z[r, c]
      frob_norm <- sqrt(sum(b_sub^2))
      if (frob_norm > lambda * w[i, j]) {
        return((1 - lambda * w[i, j] / frob_norm) * b_sub)
      } else {
        return(matrix(0, nrow = length(r), ncol = length(c)))
      }
    }
    # Initialize the final matrix
    final <- matrix(0, nrow = nrow(beta_z), ncol = ncol(beta_z))
    # Apply the function over all combinations of row and column indices
    for (i in 1:length(row_indices)) {
      for (j in 1:length(col_indices)) {
        final[row_indices[[i]], col_indices[[j]]] <- process_submatrix(row_indices[[i]], col_indices[[j]], i, j)
      }
    }
    return(final)
  }
  
  ##### full sample generation
  generate_multinomial_final<-function(J,d,p,n_sample, ind_start_X,ind_end_X,r=1/2,s=length(J),beta=NULL,R=NULL){
    obj=design_H(J,d)
    H=obj$H
    #L_sum
    ind_start_H=obj$ind_start
    ind_end_H=obj$ind_end
    #p
    L_sum=ncol(H)
    # only first order
    if (is.null(beta)) {
      beta=beta_generate(L_sum,p,ind_end=ind_end_H,s)
    }
    if (is.null(R)) {
      R = chol(cov_decay(p-1, r)) # Sigma == t(R) %*%  R
    }
    X = t(R) %*% matrix(rnorm( (p-1)*n_sample), nrow=p-1 )
    X = rbind(rep(1,n_sample), X)
    Y=multinomial_Y_generate(H, X, beta)
    return(list(X=X, Y=Y, beta=beta, R=R, H=H, ind_start_H=ind_start_H, ind_end_H=ind_end_H ))
  }
  ### no need
  generate_Poisson_final<-function(J,d,p,n_sample, ind_start_X,ind_end_X,r=1/2,s=length(J)){
    obj=design_H(J,d)
    H=obj$H
    #L_sum
    ind_start_H=obj$ind_start
    ind_end_H=obj$ind_end
    #p
    L_sum=ncol(H)
    # only first order
    beta=beta_generate(L_sum,p,ind_end=ind_end_H,s)
    R = chol(cov_decay(p-1, r)) # Sigma == t(R) %*%  R
    X = t(R) %*% matrix(rnorm( (p-1)*n_sample), nrow=p-1 )
    X = rbind(rep(1,n_sample), X)
    Y=poisson_Y_generate(H, X, beta)
    return(list(X=X, Y=Y, beta=beta, R=R, H=H, ind_start_H=ind_start_H, ind_end_H=ind_end_H ))
  }
  
  ##### Lambda max
  group_lasso_matrix <- function(beta, ind_start_H,ind_end_H,ind_start_X,ind_end_X ){
    # Create a list of row and column indices
    row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
    col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
    
    # Function to extract submatrix and compute Frobenius norm
    get_norm <- function(r, c) {
      sub_matrix <- beta[r, c, drop = FALSE]
      norm(sub_matrix, type = "F")
    }
    # Apply the function over all combinations of row and column indices
    norms <- outer(row_indices, col_indices, Vectorize(get_norm))
    return(  norms  )
  }
  
  ###  changes from here
  lambda_max_mult <- function(H,X,Y,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
    beta_zero = matrix( rep(0, ncol(H)*nrow(X) ), nrow=ncol(H) )
    # lambda_max = max( group_lasso_matrix(grad_ll_mult(beta_zero,H,X,Y),ind_start_H,ind_end_H,ind_start_X,ind_end_X ) / w )
    division_result = group_lasso_matrix(grad_ll_mult(beta_zero,H,X,Y),ind_start_H,ind_end_H,ind_start_X,ind_end_X ) / w
    # Replace Inf and -Inf values with NA
    division_result[is.infinite(division_result)] <- NA
    # Find the non-infinity largest value
    lambda_max <- max(division_result, na.rm = TRUE)
    return( lambda_max )
  }
  
  lambda_max_pois <- function(H,X,Y,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
    beta_zero = matrix( rep(0, ncol(H)*nrow(X) ), nrow=ncol(H) )
    # lambda_max = max( group_lasso_matrix(grad_ll_pois(beta_zero,H,X,Y),ind_start_H,ind_end_H,ind_start_X,ind_end_X )/w )
    division_result = group_lasso_matrix(grad_ll_pois(beta_zero,H,X,Y),ind_start_H,ind_end_H,ind_start_X,ind_end_X )/w 
    # Replace Inf and -Inf values with NA
    division_result[is.infinite(division_result)] <- NA
    # Find the non-infinity largest value
    lambda_max <- max(division_result, na.rm = TRUE)
    return(lambda_max)
  }
  
  lambda_seq<-function(lambda_max,lambda.min.ratio=0.618,nlambda=100){
    exp(seq(from = log(lambda_max),  by = log(lambda.min.ratio), length.out = nlambda ))
  }
  
  weight_generate <- function( J, d, ind_end_H, ind_end_X, weight_deep=NULL, weight_groupX=NULL, scheme_X=0 ){
    #  scheme_X = 0, w > 0
    #  scheme_X = 1, w_0,0 = 0
    #  scheme_X = 2, w_k,0=0, ||k||_0=1
    
    if ( is.null(weight_groupX) ){
      weight_groupX=rep(1,length(ind_end_X))
    }
    if ( is.null(weight_deep) ){
      weight_deep_vec=rep(1,length(ind_end_H))
    }else{
      q=length(J)
      weight_deep_vec=weight_deep[1]
      for(i in 1:d){
        weight_deep_vec=c(weight_deep_vec, rep(weight_deep[i+1], choose(q, i)) )
      }
    }
    weight_mat = outer(weight_deep_vec, weight_groupX, "*")
    if( scheme_X == 1 ){
      weight_mat[1:ind_end_H[1], 1 ] =0
    }
    if( scheme_X == 2 ){
      weight_mat[ 1:ind_end_H[ length(J)+1] , 1 ] =0
    }
    
    return(weight_mat)
  }
  
  ### changes end from here
  
  ########## Optimization functions
  ##### Backtracking Proximal Gradient Descent
  backtracking_PGD_single<-function(gradient,loss,beta_init=NULL,H,X,Y,step_size,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lambda,shrinking_gamma=0.618,epsilon=10^{-8},max_iter=1000,min_iter=10, w=1){
    #constants
    gamma=shrinking_gamma
    L_sum=ncol(H)
    #lambda
    #  nlambda=length(lam_seq)
    eta=step_size
    if(is.null(beta_init)){
      beta_loop=matrix(rep(0, ncol(H)*nrow(X) ) , nrow=ncol(H)  )
    }else{beta_loop=beta_init}
    cost=rep(0,max_iter)
    {
      z=beta_loop
      log_loss_z=loss(z, H,X,Y)
      grad=gradient(z, H,X,Y)
      repeat{
        beta_loop_new=group_lasso_proximal(z-eta*grad , w=w, lambda*eta, ind_start_H, ind_end_H, ind_start_X, ind_end_X)
        # Calculate D_loop
        D_loop = beta_loop_new - z 
        # Check step - Adjust eta if the condition is met
        log_loss_beta =loss(beta_loop_new, H, X, Y)
        #dis_loss=ll_mult(beta_loop_new, H, X, Y) - ll_mult(beta_loop, H, X, Y)
        if ( log_loss_beta <= log_loss_z + sum(grad * D_loop)+sum(D_loop^2) /(eta * 2 ) ) { break }
        eta = gamma * eta
      }
      cst_old=log_loss_beta + lambda*group_lasso(beta_loop_new ,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=w) 
      for (i in 1:max_iter) {
        z = beta_loop_new + i/(i+3)* ( beta_loop_new-beta_loop)
        beta_loop = beta_loop_new 
        grad=gradient(z, H,X,Y)
        log_loss_z=loss(z, H, X, Y)
        repeat{
          # group lasso update
          beta_loop_new=group_lasso_proximal(z-eta*grad, w=w, lambda*eta, ind_start_H, ind_end_H, ind_start_X, ind_end_X)
          # Calculate D_loop
          D_loop = beta_loop_new - z 
          # Check step - Adjust eta if the condition is met
          log_loss_beta =loss(beta_loop_new, H, X, Y)
          if ( log_loss_beta <= log_loss_z + sum(grad * D_loop)+sum(D_loop^2) /(eta * 2 ) ) { break }
          eta = gamma * eta
        }
        cost[i]=log_loss_beta+lambda*group_lasso(beta_loop_new, ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=w)  
        cst_new=cost[i]
        if (cst_old -cst_new < epsilon & sum(D_loop^2) < epsilon ){
          if (i>min_iter){
            beta_final=beta_loop_new
            break}
        }
        if (i == max_iter){beta_final=beta_loop_new}
        cst_old=cst_new
      }
      training_loss=cost
    }
    return(list(beta_final=beta_final,training_loss=training_loss[1:i], eta_end=eta))
  }
  
  backtracking_PGD_path<-function(gradient,loss,H,X,Y,step_size,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq,shrinking_gamma=0.618,epsilon=10^{-8},max_iter=100,min_iter=10,w=1){
    #constants
    gamma=shrinking_gamma
    L_sum=ncol(H)
    #lambda
    nlambda=length(lam_seq)
    #lambda_max
    ot=backtracking_PGD_single(gradient,loss,beta_init=NULL,H,X,Y,step_size,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq[1],shrinking_gamma,epsilon,max_iter,min_iter=min_iter,w)  
    eta=ot$eta_end
    beta_tensor=array( rep(0, prod(dim(ot$beta_final))*nlambda ),dim=c(ncol(H),nrow(X) ,nlambda)   )
    beta_tensor[,,1]=ot$beta_final
    #  training_loss=matrix(rep(0, max_iter*nlambda ) , nrow=max_iter  )
    lam_ind=1
    for(lambda in lam_seq[-1]){
      cost=rep(0,max_iter)
      ot=backtracking_PGD_single(gradient,loss,beta_init=beta_tensor[,,lam_ind],H,X,Y,step_size=eta/gamma,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lambda,shrinking_gamma,epsilon,max_iter,min_iter=min_iter,w)  
      eta=ot$eta_end
      lam_ind=lam_ind+1
      beta_tensor[,,lam_ind]=ot$beta_final
    }
    return( list( beta_tensor=beta_tensor ) )
  }
  
  ##### validation set tuning
  beta_validation<-function(beta_tensor,lam_seq,loss,H,X_validation,Y_validation,n_train,plot_TF=FALSE){
    nlambda=dim(beta_tensor)[3]
    #cross entropy
    cross_entropy=rep(0, nlambda)
    for(i in 1:nlambda){
      cross_entropy[i]=loss(beta_tensor[,,i],H,X_validation,Y_validation)
    }
    # Step 1: Find the minimum cross_entropy and its index
    min_cross_entropy <- min(cross_entropy)
    min_index <- which.min(cross_entropy)
    #? Step 2: Calculate the standard error of the cross_entropy values
    #ot=ll_mult_se(beta_tensor[,,min_index], H, X_validation, Y_validation)
    se <- ll_mult_se(beta_tensor[,,min_index], H, X_validation, Y_validation)$se_ll / sqrt(n_train)
    # Step 3: Find the largest lambda within one standard error of the minimum cross_entropy
    cross_entropy_threshold <- min_cross_entropy + se
    valid_lambdas <- lam_seq[ cross_entropy <= cross_entropy_threshold ]
    lambda_1se <- max(valid_lambdas)
    lambda_1se_index = min(which( cross_entropy <= cross_entropy_threshold ))
    if(plot_TF == TRUE){
      # Plotting the cross_entropy function vs lambda
      plot(log(lam_seq), cross_entropy, type = "l", xlab = "log(λ)", ylab = "Cross Entropy", main = "Cross Entropy vs λ")
      # Adding a horizontal line for min_cross_entropy + se
      abline(h = cross_entropy_threshold, col = "red", lty = 2)
      # Marking the lambda_1se, lambda.min on the plot
      points(log(lambda_1se), min(cross_entropy[lam_seq == lambda_1se]), col = "blue", pch = 3)
      #text(log(lambda_1se), 0, labels ="λ.1se" , pos = 3)
      text(log(lambda_1se), min(cross_entropy[lam_seq == lambda_1se]), labels ="λ.1se" , pos =3, col = "blue")
      points(log(lam_seq[min_index]), min_cross_entropy, col = "darkgreen", pch =4)
      text(log(lam_seq[min_index]), min_cross_entropy, labels = "λ.min " , pos = 2,col = "darkgreen")
    }
    return(list(lambda_min_index=min_index, lambda_1se_index=lambda_1se_index ))
  }
  
  backtracking_PGD_path_validation<-function(gradient,loss,H,X_train,Y_train,X_test,Y_test,step_size,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq,shrinking_gamma=0.618,epsilon=10^{-8},max_iter=100, w=1){
    #browser()
    #constants
    gamma=shrinking_gamma
    L_sum=ncol(H)
    #lambda
    nlambda=length(lam_seq)
    #lambda_max
    ot=backtracking_PGD_single(gradient,loss,beta_init=NULL,H,X_train,Y_train,step_size,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq[1],shrinking_gamma,epsilon,max_iter,w)  
    beta_tensor=array( rep(0, prod(dim(ot$beta_final))*nlambda ),dim=c(ncol(H),nrow(X_train) ,nlambda))
    beta_tensor[,,1]=ot$beta_final
    #cross entropy
    cross_entropy=rep(0, nlambda)
    cross_entropy[1]=loss(beta_tensor[,,1],H,X_test,Y_test)
    for(i in 2:nlambda){
      eta=ot$eta_end
      ot=backtracking_PGD_single(gradient,loss,beta_init=beta_tensor[,,i-1],H,X_train,Y_train,step_size=eta/gamma,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq[i],shrinking_gamma,epsilon,max_iter,w)  
      cross_entropy[i]=loss(ot$beta_final,H,X_test,Y_test)
      beta_tensor[,,i]=ot$beta_final
    }
    # Step 1: Find the minimum cross_entropy and its index
    min_cross_entropy <- min(cross_entropy)
    min_index <- which.min(cross_entropy)
    # Step 2: Calculate the standard error of the cross_entropy values
    se <- sd(cross_entropy) / sqrt(length(cross_entropy))
    # Step 3: Find the largest lambda within one standard error of the minimum cross_entropy
    cross_entropy_threshold <- min_cross_entropy + se
    valid_lambdas <- lam_seq[ cross_entropy <= cross_entropy_threshold ]
    lambda_1se <- max(valid_lambdas)
    lambda_1se_index = min(which( cross_entropy <= cross_entropy_threshold ))
    # Plotting the cross_entropy function vs lambda
    plot(log(lam_seq), cross_entropy, type = "l", xlab = "log(λ)", ylab = "Cross Entropy", main = "Cross Entropy vs λ")
    # Adding a horizontal line for min_cross_entropy + se
    abline(h = cross_entropy_threshold, col = "red", lty = 2)
    # Marking the lambda_1se, lambda.min on the plot
    points(log(lambda_1se), min(cross_entropy[lam_seq == lambda_1se]), col = "blue", pch = 3)
    #text(log(lambda_1se), 0, labels ="λ.1se" , pos = 3)
    text(log(lambda_1se), min(cross_entropy[lam_seq == lambda_1se]), labels ="λ.1se" , pos =3, col = "blue")
    points(log(lam_seq[min_index]), min_cross_entropy, col = "darkgreen", pch =4)
    text(log(lam_seq[min_index]), min_cross_entropy, labels = "λ.min " , pos = 2,col = "darkgreen")
    return(list(beta_tensor=beta_tensor,cross_entropy=cross_entropy,lambda_min_index=min_index, lambda_1se_index=lambda_1se_index ))
  }
  
  ###### tensor reshaping
  tensor_marginal <- function( J, Y, marginal_int ){
    nn=ncol(Y)
    output_mar=matrix(0, nrow = J[marginal_int], ncol=nn )
    for(i in 1:nn){
      y=Y[,i]
      y_tensor <- array(Y[,i], dim =J)
      output_mar[,i]=apply(y_tensor, marginal_int, sum)
    }
    return(output_mar)
  }
  
  ##### validation set tuning for glmnet
  beta_validation_glmnet_theta<-function(glmnet_output,X_validation,Y_validation,plot_TF=FALSE){
    loss=ll_mult_theta
    lam_seq <- glmnet_output$lambda
    lam_seq=glmnet_output$lambda
    nlambda=length(lam_seq)
    #cross entropy
    cross_entropy=rep(0, nlambda)
    for(i in 1:nlambda){
      coef_list <- coef(glmnet_output, s = lam_seq[i] )
      theta_est <- t(do.call(cbind, lapply(coef_list, as.matrix)) )
      cross_entropy[i]=loss(theta_est,X_validation,Y_validation )
    }
    # Step 1: Find the minimum cross_entropy and its index
    min_cross_entropy <- min(cross_entropy)
    min_index <- which.min(cross_entropy)
    if(plot_TF == TRUE){
      # Plotting the cross_entropy function vs lambda
      plot(log(lam_seq), cross_entropy, type = "l", xlab = "log(λ)", ylab = "Cross Entropy", main = "Cross Entropy vs λ")
      # Marking the lambda_1se, lambda.min on the plot
      points(log(lam_seq[min_index]), min_cross_entropy, col = "darkgreen", pch =4)
      text(log(lam_seq[min_index]), min_cross_entropy, labels = "λ.min " , pos = 2,col = "darkgreen")
    }
    theta.min_list <- coef(glmnet_output, s = lam_seq[min_index] )
    theta.min <- t(do.call(cbind, lapply(theta.min_list, as.matrix)) )
    return( list(lambda_min_index=min_index, theta.min=theta.min) )
  }
  
  sep_mult_beta<- function(J, X_train, Y_train, X_validation, Y_validation, plot_TF=FALSE) {
    # Get the length of J to determine the number of groups
    q <- length(J)
    # Initialize the list to store outpp$theta.min for each k
    theta_min_list <- vector("list", q)
    for (k in 1:q) {
      # Use the tensor_marginal function to separate class labels for training data
      class_Y_train_sep <- apply(tensor_marginal(J, Y_train, k), 2, which.max)
      # Use the tensor_marginal function for validation data
      class_Y_validation_sep <- tensor_marginal(J, Y_validation, k)
      # Generate the design matrix H for each group using the design_H function
      H_sep <- design_H(J[k], 1)$H
      # Fit the glmnet model and suppress potential warnings
      suppressWarnings({
        glmnet_output <- glmnet(t(X_train[-1, ]), class_Y_train_sep, family = "multinomial", type.multinomial = "grouped")
      })
      # Validate the model and get theta.min using the beta_validation_glmnet_theta function
      outpp <- beta_validation_glmnet_theta(glmnet_output, X_validation, class_Y_validation_sep, plot_TF = plot_TF)
      # Store the theta.min value in the list for each k
      theta_min_list[[k]] <- outpp$theta.min
    }
    theta_ret=matrix(0, nrow=prod(J) , ncol=nrow(X) )
    # Loop through J from 1 to the length of J using k as the loop variable
    for (k in 1:length(J)) {
      rt=prod(J[1:q<k])
      lt=prod(J[1:q>k])
      rt_vec <- matrix(rep(1, rt ), ncol = 1)
      lt_vec <- matrix(rep(1, lt ), ncol = 1)
      mid_value <- kronecker( theta_min_list[[k]], rt_vec)
      kron_result <- kronecker( lt_vec, mid_value )
      theta_ret=theta_ret+kron_result
    }
    # Return the list of theta.min values
    return(theta_ret)
  }
  
  signal_strengthen <- function(beta, delta=2, sigma=2) {
    # Create a matrix of the same size as beta
    result <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    # Apply the transformation to each element
    for (i in 1:nrow(beta)) {
      for (j in 1:ncol(beta)) {
        if (beta[i, j] > 0) {
          result[i, j] <- (beta[i, j] + delta) / sigma
        } else if (beta[i, j] < 0) {
          result[i, j] <- (beta[i, j] - delta) / sigma
        } else {
          # If beta[i, j] is zero, keep it as zero
          result[i, j] <- 0
        }
      }
    }
    return(result)
  }
  
  ######### Evaluation: Hellinger distance/Misclassification rate/F-measure and G-measure
  ##### Hellinger distance
  Hellinger_distance_ave <- function(X_test, beta, beta_est, H, theta=NULL) {
    X_mat=X_test
    # Precompute H_beta %*% X_mat and H_beta_est %*% X_mat, and apply exp
    if( is.null(theta) ){
      est_exp_mat <- exp(H %*% beta_est %*% X_mat)
    }else{
      est_exp_mat <- exp(theta %*% X_mat)
    }
    true_exp_mat <- exp(H %*% beta %*% X_mat)
    # Normalize each column by its sum to get PMFs, ensuring each column sums to 1
    true_PMFs <- sweep(true_exp_mat, 2, colSums(true_exp_mat), FUN="/")
    est_PMFs <- sweep(est_exp_mat, 2, colSums(est_exp_mat), FUN="/")
    # Calculate Hellinger distances using vectorized operations
    hellinger_dist <- sqrt(colSums((sqrt(true_PMFs) - sqrt(est_PMFs))^2)) / sqrt(2)
    # Return the average Hellinger distance
    return(mean(hellinger_dist))
  }
  
  ##### Misclassification rate
  calculate_misclassification_rate <- function(X_test, Y_test, beta_est, H, theta=NULL) {
    # Compute the exponentiated values
    if( is.null(theta) ){
      exp_mat <- exp(H %*% beta_est %*% X_test)
    }else{
      exp_mat <- exp(theta %*% X_test)
    }
    # Normalize each column to get predicted PMFs
    predicted_probabilities <- sweep(exp_mat, 2, colSums(exp_mat), FUN="/")
    # Find the index of the maximum predicted probability for each observation
    predicted_classes <- max.col( t(predicted_probabilities) )
    # Assuming Y_test is one-hot encoded, find the true class index for each observation
    true_classes <- max.col( t(Y_test) )
    # Calculate the number of misclassifications
    num_misclassified <- sum(predicted_classes != true_classes)
    # Calculate the total number of observations
    total_observations <- ncol(Y_test)
    # Calculate the misclassification rate
    misclassification_rate <- num_misclassified / total_observations
    return(misclassification_rate)
  }
  
  ##### Define a function to calculate F-measure and G-measure for variable selection
  ###   No L-Mult-theta, Sep-Mult
  F_G_Measures_old <- function(beta, beta_est, ind_start_H, ind_end_H, ind_start_X, ind_end_X) {
    # Compute group norms for the true beta coefficients using the group_lasso_matrix function
    g_norm_beta = group_lasso_matrix(beta, ind_start_H, ind_end_H, ind_start_X, ind_end_X)
    # Compute group norms for the estimated beta coefficients
    g_norm_beta_est = group_lasso_matrix(beta_est, ind_start_H, ind_end_H, ind_start_X, ind_end_X)
    # Determine the support of beta and beta_est by checking which group norms are greater than 0
    trueSupport = (g_norm_beta > 0)[-1,]
    estimatedSupport = (g_norm_beta_est > 0)[-1,]
    # Calculate the intersection of true and estimated supports by counting groups active in both
    intersection = sum(trueSupport & estimatedSupport)
    # Calculate the size of the true and estimated supports (number of active groups)
    trueSize = sum(trueSupport)
    estimatedSize = sum(estimatedSupport)
    # Calculate the F-measure, which is the harmonic mean of precision and recall
    F_measure = ifelse((trueSize + estimatedSize) > 0, (2 * intersection) / (trueSize + estimatedSize), 0)
    # Calculate the G-measure, similar to F-measure but with geometric mean
    G_measure = ifelse((trueSize * estimatedSize) > 0, (intersection) / sqrt(trueSize * estimatedSize), 0)
    # Return a list containing the F-measure and G-measure
    return(list(F_measure = F_measure, G_measure = G_measure))
  }
  
  F_G_Measures <- function(beta, beta_est, ind_start_H, ind_end_H, ind_start_X, ind_end_X) {
    # Existing code to compute group norms and supports
    g_norm_beta = group_lasso_matrix(beta, ind_start_H, ind_end_H, ind_start_X, ind_end_X)
    g_norm_beta_est = group_lasso_matrix(beta_est, ind_start_H, ind_end_H, ind_start_X, ind_end_X)
    trueSupport = (g_norm_beta > 0)[-1,]
    estimatedSupport = (g_norm_beta_est > 0)[-1,]
    
    # Calculate intersection, trueSize, and estimatedSize as before
    intersection = sum(trueSupport & estimatedSupport)
    trueSize = sum(trueSupport)
    estimatedSize = sum(estimatedSupport)
    
    # Calculate FPR and FNR
    FPR = sum(!trueSupport & estimatedSupport) / sum(!trueSupport) # FP / (FP + TN)
    FNR = sum(trueSupport & !estimatedSupport) / sum(trueSupport)  # FN / (TP + FN)
    
    # Calculate F-measure and G-measure as before
    F_measure = ifelse((trueSize + estimatedSize) > 0, (2 * intersection) / (trueSize + estimatedSize), 0)
    G_measure = ifelse((trueSize * estimatedSize) > 0, (intersection) / sqrt(trueSize * estimatedSize), 0)
    
    # Return a list containing the F-measure, G-measure, FPR, and FNR
    return(list(F_measure = F_measure, G_measure = G_measure, FPR = FPR, FNR = FNR))
  }
  
  ######### Basic functions theta version
  ll_mult_theta <- function(theta,X,Y){
    n_sa=ncol(Y)
    u_mat=theta %*% X 
    C_vec=log(colSums(exp(u_mat)))
    C=as.numeric(colSums(Y)%*%C_vec)
    return( (-sum(Y*u_mat) + C)/n_sa )
  }
  
  grad_ll_mult_theta <- function(theta,X,Y){
    expu_mat=exp(theta %*% X) 
    n_sa=ncol(Y)
    #expsum_vec= colSums(exp(u_mat))
    term_mat=(rep(1,nrow(expu_mat))%*%t(colSums( expu_mat )) )
    p_mat=expu_mat/term_mat
    dis_mat=p_mat-Y
    term=0
    for(j in 1:n_sa){
      term=term+dis_mat[,j]%*%t(X[,j])
    }
    return( term/n_sa)
  }
  
  lasso_proximal_theta <- function(Z, lambda) {
    # Apply the soft-thresholding operation element-wise to the matrix Z
    theta <- matrix(data = 0, nrow = nrow(Z), ncol = ncol(Z))  # Initialize the result matrix
    # Vectorized operation for soft-thresholding
    theta <- sign(Z) * pmax(abs(Z) - lambda, 0)
    return(theta)
  }
  
  lambda_max_lasso <- function(X, Y, w = 1) {
    # Initialize theta to zero
    theta = matrix(0, nrow = nrow(Y), ncol = nrow(X))
    # Compute the gradient of the loss function at the zero-initialized theta
    gradient_at_zero = grad_ll_mult_theta(theta, X, Y)
    lambda_max = max( sqrt(rowSums(gradient_at_zero^2))  / w)
    return(lambda_max)
  }
  
  group_lasso_proximal_theta <- function(theta_z,  lambda) {
    # Compute Frobenius norm for each row
    frob_norms <- sqrt(rowSums(theta_z^2))
    # Compute scaling factors for rows where the Frobenius norm exceeds lambda * weight
    scaling_factors <- pmax(1 - lambda / frob_norms, 0)
    # Apply scaling factor to each element of theta_z
    result <- sweep(theta_z, 1, scaling_factors, '*')
    return(result)
  }
  
  lambda_max_G_lasso <- function(X, Y, w = 1) {
    # Initialize theta to zero
    theta = matrix(0, nrow = nrow(Y), ncol = nrow(X))
    # Compute the gradient of the loss function at the zero-initialized theta
    gradient_at_zero = grad_ll_mult_theta(theta, X, Y)
    # For Lasso, lambda_max is the max absolute value of the gradient components, scaled by w
    lambda_max = max( abs(gradient_at_zero)  / w)
    return(lambda_max)
  }
  
  accelerated_proximal_gradient_descent_fixed_eta <- function(X, Y, lambda, step_size, epsilon=10^{-8},max_iter=1000,min_iter=10, theta_init=NULL,L=1) {
    if(L==1){proximal_solver=lasso_proximal_theta}else{proximal_solver=group_lasso_proximal_theta}
    eta = step_size
    if(is.null(theta_init)){
      theta = matrix(0, nrow = nrow(Y), ncol = nrow(X))  # Initialize theta
    }else{theta=theta_init}
    theta_old <- theta  # Initialize previous theta
    t_old <- 1  # Initialize momentum term
    cost=numeric(max_iter)
    for (i in 1:max_iter) {
      z = theta + ((t_old - 1) / (t_old + 2)) * (theta - theta_old)  # Calculate momentum term
      grad <- grad_ll_mult_theta(z, X, Y)  # Compute the gradient at z
      log_loss_z <- ll_mult_theta(z, X, Y)  # Compute the loss at z
      #    repeat{
      theta_new <- proximal_solver(z - eta * grad, lambda * eta)  # Update theta using proximal operator
      D_loop <- theta_new - z  # Calculate change in theta
      log_loss_theta_new <- ll_mult_theta(theta_new, X, Y)  # Compute the loss at the new theta
      if(L==1){cost[i] <- log_loss_theta_new + lambda * sum(abs(theta_new))  }else{ 
        cost[i] <- log_loss_theta_new + lambda * sum(sqrt(rowSums(theta_new^2)))  }
      # Check for convergence
      if(i>1){
        if (abs(cost[i] - cost[i-1]) < epsilon && sum(D_loop^2) < epsilon) {
          if (i > min_iter) {
            break  # Break if minimum iterations are met and convergence criteria are satisfied
            theta_old <- theta
          }
        }
      }
      # Update for the next iteration
      theta_old <- theta
      theta <- theta_new
      t_old <- t_old + 1  # Update momentum term
      cost=cost[1:i]
    }
    return(list(theta_final = theta_new, eta_end = eta,cost=cost))
  }
  
  accelerated_proximal_gradient_descent_path_fixed_eta <- function(X, Y, lam_seq, step_size, epsilon = 10^-8, max_iter = 100, min_iter = 10, L=1) {
    #browser()  # Start debugging here
    # Determine dimensions for the tensor
    nlambda <- length(lam_seq)
    nfeatures <- nrow(Y)  # Assuming theta has the same number of rows as X has features
    nresponses <- nrow(X)  # Assuming theta has the same number of columns as Y has responses
    # Initialize the 3D tensor to store beta matrices for each lambda
    theta_tensor <- array(0, dim = c(nfeatures, nresponses, nlambda))
    # Initialize step size and other parameters
    eta <- step_size
    theta_init <- NULL  # Initialize theta for the first lambda
    for (lambda_index in 1:nlambda) {
      lambda <- lam_seq[lambda_index]
      # Run accelerated proximal gradient descent for the current lambda
      optimization_result <- accelerated_proximal_gradient_descent_fixed_eta(X, Y, lambda, eta, epsilon, max_iter, min_iter, theta_init,L=L)
      # Store the final theta in the tensor
      theta_tensor[,,lambda_index] <- optimization_result$theta_final
      # Use the final theta from this lambda as the initial theta for the next lambda
      theta_init <- optimization_result$theta_final
      # Update eta for the next lambda based on the final eta from the current optimization
      # eta <- optimization_result$eta_end/gamma
    }
    return(theta_tensor)
  }
  
  beta_validation_theta<-function(theta_tensor,lam_seq,X_validation,Y_validation,plot_FT=FALSE){
    nlambda=dim(theta_tensor)[3]
    #cross entropy
    cross_entropy=rep(0, nlambda)
    for(i in 1:nlambda){
      cross_entropy[i]=ll_mult_theta(theta_tensor[,,i],X_validation,Y_validation)
    }
    # Step 1: Find the minimum cross_entropy and its index
    min_cross_entropy <- min(cross_entropy)
    min_index <- which.min(cross_entropy)
    if(plot_FT == TRUE){
      # Plotting the cross_entropy function vs lambda
      plot(log(lam_seq), cross_entropy, type = "l", xlab = "log(λ)", ylab = "Cross Entropy", main = "Cross Entropy vs λ")
      # Marking the lambda_1se, lambda.min on the plot
      points(log(lam_seq[min_index]), min_cross_entropy, col = "darkgreen", pch =4)
      text(log(lam_seq[min_index]), min_cross_entropy, labels = "λ.min " , pos = 2,col = "darkgreen")
    }
    return(lambda_min_index=min_index )
  }
  
  ######### Basic functions \theta version
  
  ll_mult_theta <- function(theta,X,Y){
    n_sa=ncol(Y)
    u_mat=theta %*% X 
    C_vec=log(colSums(exp(u_mat)))
    C=as.numeric(colSums(Y)%*%C_vec)
    return( (-sum(Y*u_mat) + C)/n_sa )
  }
  
  grad_ll_mult_theta <- function(theta,X,Y){
    expu_mat=exp(theta %*% X) 
    n_sa=ncol(Y)
    #expsum_vec= colSums(exp(u_mat))
    term_mat=(rep(1,nrow(expu_mat))%*%t(colSums( expu_mat )) )
    p_mat=expu_mat/term_mat
    dis_mat=p_mat-Y
    term=0
    for(j in 1:n_sa){
      term=term+dis_mat[,j]%*%t(X[,j])
    }
    return( term/n_sa)
  }
  
  lasso_proximal_theta <- function(Z, lambda) {
    # Apply the soft-thresholding operation element-wise to the matrix Z
    theta <- matrix(data = 0, nrow = nrow(Z), ncol = ncol(Z))  # Initialize the result matrix
    # Vectorized operation for soft-thresholding
    theta <- sign(Z) * pmax(abs(Z) - lambda, 0)
    return(theta)
  }
  
  lambda_max_lasso <- function(X, Y, w = 1) {
    # Initialize theta to zero
    theta = matrix(0, nrow = nrow(Y), ncol = nrow(X))
    # Compute the gradient of the loss function at the zero-initialized theta
    gradient_at_zero = grad_ll_mult_theta(theta, X, Y)
    lambda_max = max( sqrt(rowSums(gradient_at_zero^2))  / w)
    return(lambda_max)
  }
  
  group_lasso_proximal_theta <- function(theta_z,  lambda) {
    # Compute Frobenius norm for each row
    frob_norms <- sqrt(rowSums(theta_z^2))
    # Compute scaling factors for rows where the Frobenius norm exceeds lambda * weight
    scaling_factors <- pmax(1 - lambda / frob_norms, 0)
    # Apply scaling factor to each element of theta_z
    result <- sweep(theta_z, 1, scaling_factors, '*')
    return(result)
  }
  
  lambda_max_G_lasso <- function(X, Y, w = 1) {
    # Initialize theta to zero
    theta = matrix(0, nrow = nrow(Y), ncol = nrow(X))
    # Compute the gradient of the loss function at the zero-initialized theta
    gradient_at_zero = grad_ll_mult_theta(theta, X, Y)
    # For Lasso, lambda_max is the max absolute value of the gradient components, scaled by w
    lambda_max = max( abs(gradient_at_zero)  / w)
    return(lambda_max)
  }
  
  accelerated_proximal_gradient_descent_fixed_eta <- function(X, Y, lambda, step_size, epsilon=10^{-8},max_iter=1000,min_iter=10, theta_init=NULL,L=1) {
    if(L==1){proximal_solver=lasso_proximal_theta}else{proximal_solver=group_lasso_proximal_theta}
    eta = step_size
    if(is.null(theta_init)){
      theta = matrix(0, nrow = nrow(Y), ncol = nrow(X))  # Initialize theta
    }else{theta=theta_init}
    theta_old <- theta  # Initialize previous theta
    t_old <- 1  # Initialize momentum term
    cost=numeric(max_iter)
    for (i in 1:max_iter) {
      z = theta + ((t_old - 1) / (t_old + 2)) * (theta - theta_old)  # Calculate momentum term
      grad <- grad_ll_mult_theta(z, X, Y)  # Compute the gradient at z
      log_loss_z <- ll_mult_theta(z, X, Y)  # Compute the loss at z
      #    repeat{
      theta_new <- proximal_solver(z - eta * grad, lambda * eta)  # Update theta using proximal operator
      D_loop <- theta_new - z  # Calculate change in theta
      log_loss_theta_new <- ll_mult_theta(theta_new, X, Y)  # Compute the loss at the new theta
      if(L==1){cost[i] <- log_loss_theta_new + lambda * sum(abs(theta_new))  }else{ 
        cost[i] <- log_loss_theta_new + lambda * sum(sqrt(rowSums(theta_new^2)))  }
      # Check for convergence
      if(i>1){
        if (abs(cost[i] - cost[i-1]) < epsilon && sum(D_loop^2) < epsilon) {
          if (i > min_iter) {
            break  # Break if minimum iterations are met and convergence criteria are satisfied
            theta_old <- theta
          }
        }
      }
      # Update for the next iteration
      theta_old <- theta
      theta <- theta_new
      t_old <- t_old + 1  # Update momentum term
      cost=cost[1:i]
    }
    return(list(theta_final = theta_new, eta_end = eta,cost=cost))
  }
  
  accelerated_proximal_gradient_descent_path_fixed_eta <- function(X, Y, lam_seq, step_size, epsilon = 10^-8, max_iter = 100, min_iter = 10, L=1) {
    #browser()  # Start debugging here
    # Determine dimensions for the tensor
    nlambda <- length(lam_seq)
    nfeatures <- nrow(Y)  # Assuming theta has the same number of rows as X has features
    nresponses <- nrow(X)  # Assuming theta has the same number of columns as Y has responses
    # Initialize the 3D tensor to store beta matrices for each lambda
    theta_tensor <- array(0, dim = c(nfeatures, nresponses, nlambda))
    # Initialize step size and other parameters
    eta <- step_size
    theta_init <- NULL  # Initialize theta for the first lambda
    for (lambda_index in 1:nlambda) {
      lambda <- lam_seq[lambda_index]
      # Run accelerated proximal gradient descent for the current lambda
      optimization_result <- accelerated_proximal_gradient_descent_fixed_eta(X, Y, lambda, eta, epsilon, max_iter, min_iter, theta_init,L=L)
      # Store the final theta in the tensor
      theta_tensor[,,lambda_index] <- optimization_result$theta_final
      # Use the final theta from this lambda as the initial theta for the next lambda
      theta_init <- optimization_result$theta_final
      # Update eta for the next lambda based on the final eta from the current optimization
      # eta <- optimization_result$eta_end/gamma
    }
    return(theta_tensor)
  }
  
  beta_validation_theta<-function(theta_tensor,lam_seq,X_validation,Y_validation,plot_FT=FALSE){
    nlambda=dim(theta_tensor)[3]
    #cross entropy
    cross_entropy=rep(0, nlambda)
    for(i in 1:nlambda){
      cross_entropy[i]=ll_mult_theta(theta_tensor[,,i],X_validation,Y_validation)
    }
    # Step 1: Find the minimum cross_entropy and its index
    min_cross_entropy <- min(cross_entropy)
    min_index <- which.min(cross_entropy)
    if(plot_FT == TRUE){
      # Plotting the cross_entropy function vs lambda
      plot(log(lam_seq), cross_entropy, type = "l", xlab = "log(λ)", ylab = "Cross Entropy", main = "Cross Entropy vs λ")
      # Marking the lambda_1se, lambda.min on the plot
      points(log(lam_seq[min_index]), min_cross_entropy, col = "darkgreen", pch =4)
      text(log(lam_seq[min_index]), min_cross_entropy, labels = "λ.min " , pos = 2,col = "darkgreen")
    }
    return(lambda_min_index=min_index )
  }
  
}

### loop for all different problem settings
for(ind_loop in 1:nrow(params)  ){
  p <- params[ind_loop,1]
  Model <- params[ind_loop,2]
  n <- params[ind_loop,3]
  ### generate initial beta
  {
    scheme_X=2
    scheme=Model
    n_train=n
    # p is defined above
    # initial seed for parameter generation
    set.seed(7777)
    J=c(2,2,2,3)
    d=4
    
    n_validation=1000
    n_test=10000
    #n_sample=n_train+n_validation+n_test
    ### generation
    ind_end_X = 1:p
    ind_start_X = c(1, 1+ind_end_X[-length(ind_end_X)]) 
    
    ### round 1 for R, L_sum
    output_round1=generate_multinomial_final(J,d,p,2,ind_start_X,ind_end_X,r=1/2,s=length(J),beta=NULL,R=NULL)
    R=output_round1$R
    H=output_round1$H
    L_sum=ncol(H)
    ind_start_H <- output_round1$ind_start_H
    ind_end_H <- output_round1$ind_end_H
    
    ind_start_X_FG=1:p
    ind_end_X_FG=1:p
    
    ### round 2 for scheme beta generation
    set.seed(7777)
    beta=beta_generate_scheme(J,d,L_sum,p,ind_start_H,ind_end_H,unif_range=c(-2,2),scheme=scheme)
    beta=signal_strengthen(beta, delta=2, sigma = 2)
    lambda.min.ratio_seq = c(0.94, 0.94, 0.94, 0.94, 0.94)
  }
  ### 
  
  ### loop ID for MC: nreps/num_parallel per batch
  for( i in ID_index+(0:(nreps/num_parallel-1) )*num_parallel ){
    
    {
      ii=1
      nlam=200
      H_dis=matrix(0, nrow=1, ncol=6)
      mis_class=matrix(0, nrow=1, ncol=7)
      F_measure=matrix(0, nrow=1, ncol=4)
      G_measure=matrix(0, nrow=1, ncol=4)
      FPR_mat=matrix(0, nrow=1, ncol=4)
      FNR_mat=matrix(0, nrow=1, ncol=4)
     
        ### seed in training-validation-testing
        set.seed(i)
        # Scheme related Multinomial model generation
        outpt_mult <- generate_multinomial_final(J, d, p, n_train, ind_start_X, ind_end_X, r = 1/2, s = length(J), beta = beta, R = R)
        X_train <- outpt_mult$X
        Y_train <- outpt_mult$Y
        # L_sum = ncol(H)
        epsilon <- 10^{-4}
        outpt_mult <- generate_multinomial_final(J, d, p, n_validation+n_test, ind_start_X, ind_end_X, r = 1/2, s = length(J), beta = beta, R = R)
        X=outpt_mult$X
        Y=outpt_mult$Y
        X_validation <- X[, 1:(n_validation  )]
        Y_validation <- Y[, 1:(n_validation  )]
        X_test <- X[, (1 + n_validation  ):(n_validation+n_test)]
        Y_test <- Y[, (1 + n_validation  ):(n_validation+n_test)]
        
        ############# Lasso
        ind_end_X = 1:p
        ind_start_X = c(1, 1+ind_end_X[-length(ind_end_X)]) 
        {weight_deep_log_ratio=0
          weight_deep=exp(seq(from=0, by=weight_deep_log_ratio, length.out=d+1 ) )
          weight_groupX_x_w=1
          weight_groupX=c(1, rep(weight_groupX_x_w,length(ind_end_X)-1 ) )
          weight = weight_generate( J, d, ind_end_H, ind_end_X, weight_deep, weight_groupX, scheme_X=scheme_X ) }
        ### 1. L-Mult
        lambda.min.ratio=lambda.min.ratio_seq[1]
        lambda_mult_max=lambda_max_mult(H,X_train,Y_train,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=weight )
        lam_seq_mult=lambda_seq(lambda_mult_max,lambda.min.ratio=lambda.min.ratio, nlambda = nlam)
        eta_mult=prod(J)*n_train/norm(X_train,"2")^2 
        # training
        beta_path_L_mult=backtracking_PGD_path(grad_ll_mult,ll_mult,H,X=X_train,Y=Y_train,step_size=eta_mult,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq_mult,shrinking_gamma=0.618,epsilon=epsilon,max_iter=1000, w=weight)
        #L_Mult_tensor=beta_path_L_mult$beta_tensor
        beta_tensor=beta_path_L_mult$beta_tensor
        # validation
        beta_val=beta_validation(beta_tensor,lam_seq_mult,ll_mult,H,X_validation,Y_validation,n_train=n_train,plot_TF=FALSE)
        # testing set evaluation
        H_dis[ii,1]=Hellinger_distance_ave(X_test, beta, beta_tensor[,,beta_val$lambda_min_index], H) 
        mis_class[ii,1]=calculate_misclassification_rate(X_test, Y_test,  beta_tensor[,,beta_val$lambda_min_index], H)
        FG=F_G_Measures(beta, beta_tensor[,,beta_val$lambda_min_index], ind_start_H, ind_end_H, ind_start_X_FG, ind_end_X_FG)
        F_measure[ii,1]=FG$F_measure
        G_measure[ii,1]=FG$G_measure
        FPR_mat[ii,1]=FG$FPR
        FNR_mat[ii,1]=FG$FNR
        
        ### 2. L-Pois
        lambda.min.ratio=lambda.min.ratio_seq[2]
        lambda_pois_max=lambda_max_pois(H,X_train,Y_train,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=weight)
        lam_seq_pois=lambda_seq(lambda_pois_max,lambda.min.ratio=lambda.min.ratio,nlambda=nlam)
        eta_pois=n_train/norm(X_train,"2")^2 
        # training
        beta_path_L_pois=backtracking_PGD_path(grad_ll_pois,ll_pois,H,X=X_train,Y=Y_train,step_size=eta_pois,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq_pois,shrinking_gamma=0.618,epsilon=epsilon,max_iter=1000, w=weight)
        beta_tensor=beta_path_L_pois$beta_tensor
        # validation
        beta_val=beta_validation(beta_tensor,lam_seq_pois,ll_mult,H,X_validation,Y_validation,n_train=n_train,plot_TF=FALSE)
        # testing set evaluation
        H_dis[ii,2]=Hellinger_distance_ave(X_test, beta, beta_tensor[,,beta_val$lambda_min_index], H) 
        mis_class[ii,2]=calculate_misclassification_rate(X_test, Y_test,  beta_tensor[,,beta_val$lambda_min_index], H)
        FG=F_G_Measures(beta, beta_tensor[,,beta_val$lambda_min_index], ind_start_H, ind_end_H, ind_start_X_FG, ind_end_X_FG)
        F_measure[ii,2]=FG$F_measure
        G_measure[ii,2]=FG$G_measure
        FPR_mat[ii,2]=FG$FPR
        FNR_mat[ii,2]=FG$FNR
        
        ############# group
        ind_end_X = c(1, p)
        ind_start_X = c(1, 2)
        {weight_deep_log_ratio=0
          weight_deep=exp(seq(from=0, by=weight_deep_log_ratio, length.out=d+1 ))
          weight_groupX_x_w=1
          weight_groupX=c(1, rep(weight_groupX_x_w,length(ind_end_X)-1 ) )
          weight = weight_generate( J, d, ind_end_H, ind_end_X, weight_deep, weight_groupX, scheme_X=scheme_X ) }
        
        ### 3. G-Mult
        lambda.min.ratio=lambda.min.ratio_seq[3]
        lambda_mult_max=lambda_max_mult(H,X_train,Y_train,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w = weight )
        # training
        lam_seq_mult=lambda_seq(lambda_mult_max,lambda.min.ratio=lambda.min.ratio,nlambda=nlam)
        eta_mult=prod(J)*n_train/norm(X_train,"2")^2 
        # validation
        beta_path_G_mult=backtracking_PGD_path(grad_ll_mult,ll_mult,H,X=X_train,Y=Y_train,step_size=eta_mult,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq_mult,shrinking_gamma=0.618,epsilon=epsilon,max_iter=1000, w = weight )
        beta_tensor=beta_path_G_mult$beta_tensor
        # testing set evaluation
        beta_val=beta_validation(beta_tensor,lam_seq_mult,ll_mult,H,X_validation,Y_validation,n_train=n_train,plot_TF=FALSE)
        H_dis[ii,3]=Hellinger_distance_ave(X_test, beta, beta_tensor[,,beta_val$lambda_min_index], H) 
        mis_class[ii,3]=calculate_misclassification_rate(X_test, Y_test,  beta_tensor[,,beta_val$lambda_min_index], H)
        FG=F_G_Measures(beta, beta_tensor[,,beta_val$lambda_min_index], ind_start_H, ind_end_H, ind_start_X_FG, ind_end_X_FG)
        F_measure[ii,3]=FG$F_measure
        G_measure[ii,3]=FG$G_measure
        FPR_mat[ii,3]=FG$FPR
        FNR_mat[ii,3]=FG$FNR
        
        ### 4. G-Pois
        lambda.min.ratio=lambda.min.ratio_seq[4]
        lambda_pois_max=lambda_max_pois(H,X_train,Y_train,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w = weight )
        lam_seq_pois=lambda_seq(lambda_pois_max,lambda.min.ratio=lambda.min.ratio,nlambda=nlam)
        eta_pois=n_train/norm(X_train,"2")^2 
        # training
        beta_path_G_pois=backtracking_PGD_path(grad_ll_pois,ll_pois,H,X=X_train,Y=Y_train,step_size=eta_pois,ind_start_H,ind_end_H,ind_start_X,ind_end_X,lam_seq_pois,shrinking_gamma=0.618,epsilon=epsilon,max_iter=1000, w = weight )
        beta_tensor=beta_path_G_pois$beta_tensor
        # validation
        beta_val=beta_validation(beta_tensor,lam_seq_pois,ll_mult,H,X_validation,Y_validation,n_train=n_train,plot_TF=FALSE)
        # testing set evaluation
        H_dis[ii,4]=Hellinger_distance_ave(X_test, beta, beta_tensor[,,beta_val$lambda_min_index], H) 
        mis_class[ii,4]=calculate_misclassification_rate(X_test, Y_test,  beta_tensor[,,beta_val$lambda_min_index], H)
        FG=F_G_Measures(beta, beta_tensor[,,beta_val$lambda_min_index], ind_start_H, ind_end_H, ind_start_X_FG, ind_end_X_FG)
        F_measure[ii,4]=FG$F_measure
        G_measure[ii,4]=FG$G_measure
        FPR_mat[ii,4]=FG$FPR
        FNR_mat[ii,4]=FG$FNR
        
        ### 5 L_Mult_theta: Group-penalized multinomial log-linear model
        lamb_max=lambda_max_G_lasso(X_train, Y_train, w = 1)
        lam_seq = lambda_seq(lamb_max,lambda.min.ratio=lambda.min.ratio_seq[5],nlambda=nlam)
        eta_theta=2*n_train/norm(X_train,"2")^2 
        theta_tensor_G=accelerated_proximal_gradient_descent_path_fixed_eta(X_train, Y_train, lam_seq, step_size=eta_theta, epsilon =epsilon, max_iter = 1000, min_iter=5,L=2)
        theta.min_ind_G=beta_validation_theta(theta_tensor_G,lam_seq,X_validation,Y_validation,plot_FT=FALSE)
        # testing set evaluation
        H_dis[ii,5]=Hellinger_distance_ave(X_test, beta, beta_est=beta, H, theta=theta_tensor_G[,,theta.min_ind_G] )
        mis_class[ii,5]=calculate_misclassification_rate(X_test, Y_test, beta_est=beta_est, H, theta=theta_tensor_G[,,theta.min_ind_G] ) 
        ###F/G_Measures  no need
        ### 6. Sep-Mult
        ott=sep_mult_beta(J, X_train, Y_train, X_validation, Y_validation, plot_TF=FALSE)
        H_dis[ii,6]=Hellinger_distance_ave(X_test, beta, beta_est=beta, H, theta=ott)
        mis_class[ii,6]=calculate_misclassification_rate(X_test, Y_test, beta_est=beta_est, H, theta=ott) 
        
        ### 7. Oracle
        # testing set evaluation
        mis_class[ii,7]=calculate_misclassification_rate(X_test, Y_test,  beta, H)
    }
    
    # -----------------------------------------------
    # Set location and filename for saving results
    # -----------------------------------------------
    savename <- paste("/home/shenx/zhao1118/AD-simulation/structure_learning_2/data_CE/Model", Model, "_p", p, "_n", n, "_repID", i, ".RDS", sep="")

    ### result list
    #will save this as RDS file later
    Results <- list(
      Hellinger_dist = H_dis,
      Misclassification_rate = mis_class,
      F_measure = F_measure,
      G_measure = G_measure, 
      FPR_mat = FPR_mat,
      FNR_mat = FNR_mat,
      Model=Model,
      p=p,
      n=n,
      rep_ID=i
    )

    # ---------------------------------------------
    # save results
    # ---------------------------------------------
    saveRDS(Results, file=savename)
  }
  
  
}




