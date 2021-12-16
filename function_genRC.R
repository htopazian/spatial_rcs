
#############################
# Algorithm to generate estimates of reproduction numbers, unobserved cases (epsilon edges), 
#' and second order AIC estimates for a particular timeseries, distance matrix and set of parameters
#'
#' @param tmat time series matrix
#' @param dmat distance matrix
#' @param DataType distance data type - euclidian, friction (malaria atlas project travel times), binary (whether two things are similar or different, 1/0)
#' @param SpatialKernel nil (no spatial information, just time), exponential (with space), gaussian
#' @param fixed vector of fixed parameters - alpha, delta, epsilon
#' @param alpha prior mean and variance for serial interval parameter alpha (temporal), if fixed the fixed value
#' @param delta prior mean and variance for distance function parameter delta (spatial), if fixed the fixed value
#' @param epsilon prior mean and variance for epsilon edge (assumption that each case has an unobserved source of infection), if fixed the fixed value. Higher values = higher likelihood of unobserved sources of infection. 
#'
#' @return
#' @export
#'
#' @examples


spatialnetrate <-
  function(tp = tmat,
           dp = dmat,
           DataType = euclidian,
           SpatialKernel = "nil",
           fixed = "epsilon",
           alpha = c(0.003, 0.01),
           delta = c(0.01, 0.01),
           epsilon = 1e-3) {
    
    library(tensorflow)
    library(reticulate)
    tfp <- reticulate::import("tensorflow_probability", convert = FALSE)
    
    
    # Data defining section --------------------------
    tf$reset_default_graph() 
    
    N = ncol(tmat) # number of infectors
    Q = nrow(tmat) # number of infectees
    
    tp = tf$placeholder(shape = shape(Q, N), dtype = tf$float32) # time matrix
    dp = tf$placeholder(shape = shape(Q, N), dtype = tf$float32) # distance matrix
    
    # A defining section
    if (fixed == "alpha") {
      A <- alpha
      A <- tf$constant(alpha, shape = shape(Q, N), dtype = tf$float32) # set as a constant
      
    } else{ # estimate if alpha is not fixed
      A <-
        tf$nn$relu(tf$Variable(  # define characteristics of what alpha could be
          tf$constant(alpha[1], shape = shape(Q, N), dtype = tf$float32),
          name = 'alpha'
        )) # length scale for space
      A <- # define a minimum alpha value = 0.01
        tf$minimum(A, tf$constant(0.01, shape = shape(Q, N), dtype = tf$float32))
      A <- A * tf$sign(tp)
      Ad <- # distribution of A = the prior that you are specifying
        tfp$distributions$TruncatedNormal(
          loc = alpha[1], # mean
          scale = alpha[2], # scale
          high = 0.01, # upper bounds
          low = 0.0 # lower bounds
        ) 
      
    } 
# each pairwise relationship can have a different alpha = a diff serial interval within specified bounds
    
    # D defining section
    if (fixed == "delta") {
      D <- tf$nn$relu(tf$constant(delta[1], shape = shape(1L), dtype = tf$float32))
      
    } else{
      D <-
        tf$nn$relu(tf$Variable(
          tf$constant(delta[1], shape = shape(1L), dtype = tf$float32),
          name = 'distance_scale'
        )) # length scale for space
      Dd <-
        tfp$distributions$TruncatedNormal(
          loc = delta[1],
          scale = delta[2],
          low = 0.00001,
          high = 0.5
        ) 
      
      
    } # just estimating one value here, not a matrix like A
    
    
    # E defining section
    if (fixed == "epsilon") {
      E <- tf$nn$relu(tf$Variable(tf$constant(epsilon,shape=shape(Q),dtype=tf$float32),name='epsilon_edge'))
      E <- tf$minimum(E,0.2) # has to be below or equal to 0.02
      
    } else{
      E <-
        tf$nn$relu(tf$Variable(
          tf$constant(epsilon[1], shape = shape(Q), dtype = tf$float32),
          name = 'epsilon_edge'
        )) # length scale for space
      Ed <-
        tfp$distributions$TruncatedNormal(epsilon[1],
                                          scale = epsilon[2],
                                          high = 1, # bound
                                          low = 0.000000000001)
    }
    if(fixed == "nil"){
      print("no fixed parameters - identifiability could be an issue")
    } # not a matrix - each non-imported case can have a different value, for each case, ask what is the likelihood that there was an unobserved infector? 
    
    
    # Spatial kernel section --------------------------
    # hazard and survival
    if (SpatialKernel == "exponential") {
      H <- A * tp * tf$exp(-dp * D)*D
      S <- tf$negative(0.5 * A * tp * tp) - tf$log(D)
      
    } else if (SpatialKernel == "gaussian"){
      H <- (2*tf$sqrt(D)*A*tp*tf$exp(-D*(dp*dp)))/(tf$sqrt(pi))
      S <- tf$negative(0.5 * A * tp * tp) + tf$log((tf$sqrt(pi))/(2*tf$sqrt(D)))
      
    } else { # 'nil'
      H <- A * tp
      S<-tf$negative(0.5 * A * tp * tp)
      
    }
    
    # Tensorflow code --------------------------
    # negative log likelihoods
    nll1 <- # hazard
      tf$reduce_sum(log(tf$reduce_sum(H, 1L) + E)) # first loop and second loop
  
    nll2 <- # survival
      tf$reduce_sum(tf$reduce_sum(S, 1L) - E) # first loop and second loop
    
    # priors
    if(fixed != "epsilon"){
      log_prior <-
        tf$negative(tf$reduce_sum(Ad$log_prob(A)) + tf$reduce_sum(Ed$log_prob(E)))
      
    } else{ 
      log_prior <-
        tf$negative(tf$reduce_sum(Ad$log_prob(A))) 
   
    }
    
    nll <- tf$negative(nll1 + nll2)
    # posterior
    log_posterior <- tf$add(nll, log_prior)
    
    F <- H * tf$exp(S) # transmission likelihood is hazard * survival
    
    
    # BFGS optimising --------------------------
    # minimizing posterior log likelihood
    optimizer = tf$contrib$opt$ScipyOptimizerInterface(
      log_posterior,
      method = 'L-BFGS-B',
      options = dict('maxiter' = 2000000L,
                     'disp' = TRUE)
    )
    sess = tf$Session()
    init = tf$global_variables_initializer()
    sess$run(init) # initialise session

    optimizer$minimize(session = sess,
                       feed_dict = dict(tp = tmat,
                                        dp = dmat))

    Ds<-sess$run(D)
    value <- sess$run(list(nll, log_prior, log_posterior),
                      feed_dict = dict(tp = tmat,
                                       dp = dmat))
    print(paste(value[[1]], value[[2]], value[[3]]))
   
     num<-nrow(tmat)
    
    if(SpatialKernel == "nil"){
      K<-num+1
    }else{
      K <- num+2
    }
    
    AIC <- 2*value[[1]] + 2*K +(2*K*(K+1)/(num-K-1))
    
    # alpha parameter
    As <- sess$run(list(A), feed_dict = dict(tp = tmat,
                                             dp = dmat))[[1]]
    meanA <- mean(As[As != 0])
    par(mfrow = c(1, 4))
    Fs <- sess$run(list(F), feed_dict = dict(tp = tmat,
                                             dp = dmat))[[1]]
    Ds <- sess$run(D)
    
    if(fixed == "epsilon"){
      eps_edge <- sess$run(E)
    } else{
      eps_edge <- sess$run(E)
    }
    
    # likelihood
    Fs <- cbind(Fs, eps_edge)
    
    Fs = sweep(Fs, 1, rowSums(Fs), FUN = '/')
    
    # R - for a particular case, look at all possible infections arising from the case and then sum across the matrix. likelihood = 0 when symptom onset occurs before candidate infector could transmit
    Rt <- colSums(Fs, na.rm=TRUE)[1:(ncol(Fs) - 1)]
    
    # plotting space-time likelihood and kernels
    x=seq(15,165,length.out=100)
    d=seq(1,200,length.out=100)
    z=expand.grid(x,d)
    r= as.numeric(Ds)
    
    return(list(Rt, Ds, As, eps_edge, AIC))
    
  }

