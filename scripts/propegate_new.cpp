// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>
#include <gsl/gsl_rng.h>
#include <gsl_randist.h>
#include <unistd.h>            // getpid


using namespace Rcpp;


#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK

// [[Rcpp::export]]
arma::mat v(arma::colvec a) {
  return a*a.t();}


// [[Rcpp::export]]
List filter_AF_pair_C(NumericVector A_t1_s1, NumericVector A_t2_s1, NumericVector A_t1_s2, NumericVector A_t2_s2){
  LogicalVector filt1 = (A_t1_s1>0.05 & A_t1_s1<0.95);
  LogicalVector filt2 = (A_t1_s2>0.05 & A_t1_s2<0.95);
  LogicalVector filt = filt1&filt2;

  
  
  return(List::create(A_t1_s1[filt], A_t2_s1[filt], A_t1_s2[filt], A_t2_s2[filt]));
}


// [[Rcpp::export]]
NumericVector freq_to_fc_C(NumericVector A_t1, NumericVector A_t2, NumericVector a_t1, NumericVector a_t2){
  return( ( ( pow((A_t1 - A_t2),2) / ((A_t1 + A_t2)/2 - A_t1 * A_t2) ) + ( pow((a_t1-a_t2),2) / ((a_t1+a_t2)/2 - a_t1*a_t2) ) ) / 2 );
}

// [[Rcpp::export]]
NumericVector freq_to_fc_sqrt_C(NumericVector A_t1, NumericVector A_t2, NumericVector a_t1, NumericVector a_t2){
  return( ( ( (A_t1 - A_t2) / sqrt(((A_t1 + A_t2)/2 - A_t1 * A_t2) )) )) ;
}

// [[Rcpp::export]]
List test_fn(arma::cube pop){
  List out(1);
  out[0] = pop;
  return out;
}

arma::rowvec rmn(unsigned int N, arma::rowvec p, gsl_rng* r){
    int K = p.n_elem;
    //arma::rowvec x(K, arma::fill::zeros);
    int x[K];
    gsl_ran_multinomial(r, K, N, p.begin(), (unsigned int *) x);
    // Rcout << p;
    // Rcout <<"\n";
    // Rcout << N;
    // Rcout <<"\n";
    // Rcout << *x;
    // Rcout <<"\n";
    // Rcout <<"\n";

    arma::rowvec o(K, arma::fill::zeros);
    for(int ii=0; ii<K; ++ii ){
        o[ii] = *(x+ii);
    }
    return o;             // return results vector
}


// [[Rcpp::export]]
arma::mat gsl_mmm(int S,  arma::mat P){
    //  << P;
    int j;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    long seed = rand()/(((double)RAND_MAX + 1)/10000000) * getpid();
    gsl_rng_set (r, seed);
    arma::mat X = arma::zeros<arma::mat>(P.n_rows, P.n_cols);
    for(j=0; j<X.n_rows; j++){
        X.row(j) = rmn(S, P.row(j), r);
    }
    gsl_rng_free (r);
    return X;
}

// [[Rcpp::export]]
Rcpp::List propegate_population_C(List pop, arma::vec N_v, arma::mat M, double c, int t_add){


  int loci = as<int>(pop["loci"]);
  int last_pre = as<int>(pop["t_curr"]);
  int n_subpops = as<int>(pop["n_subpops"]);

  arma::cube pAB_pre = as<arma::cube>(pop["pAB_pre"]);
  arma::cube pAb_pre = as<arma::cube>(pop["pAb_pre"]);
  arma::cube paB_pre = as<arma::cube>(pop["paB_pre"]);
  arma::cube pab_pre = as<arma::cube>(pop["pab_pre"]);

  arma::cube pAB_post = as<arma::cube>(pop["pAB_post"]);
  arma::cube pAb_post = as<arma::cube>(pop["pAb_post"]);
  arma::cube paB_post = as<arma::cube>(pop["paB_post"]);
  arma::cube pab_post = as<arma::cube>(pop["pab_post"]);

  arma::cube pA_pre = as<arma::cube>(pop["pA_pre"]);
  arma::cube pA_post = as<arma::cube>(pop["pA_post"]);

  arma::cube pB_pre = as<arma::cube>(pop["pB_pre"]);
  arma::cube pB_post = as<arma::cube>(pop["pB_post"]);

  arma::cube D_pre = as<arma::cube>(pop["D_pre"]);
  arma::cube D_post = as<arma::cube>(pop["D_post"]);
  arma::cube r_pre = as<arma::cube>(pop["r_pre"]);
  arma::cube r_post = as<arma::cube>(pop["r_post"]);
  
  

  //


  pAB_pre.resize( loci, n_subpops, last_pre + t_add );
  pAb_pre.resize( loci, n_subpops, last_pre + t_add );
  paB_pre.resize( loci, n_subpops, last_pre + t_add );
  pab_pre.resize( loci, n_subpops, last_pre + t_add );

  pAB_post.resize( loci, n_subpops, last_pre + t_add -1);
  pAb_post.resize( loci, n_subpops, last_pre + t_add -1);
  paB_post.resize( loci, n_subpops, last_pre + t_add -1);
  pab_post.resize( loci, n_subpops, last_pre + t_add -1);

  pA_pre.resize( loci, n_subpops, last_pre + t_add );
  pA_post.resize( loci, n_subpops, last_pre + t_add -1);

  pB_pre.resize( loci, n_subpops, last_pre + t_add );
  pB_post.resize( loci, n_subpops, last_pre + t_add -1);

  D_pre.resize( loci, n_subpops, last_pre + t_add );
  D_post.resize( loci, n_subpops, last_pre + t_add -1);
  r_pre.resize( loci, n_subpops, last_pre + t_add );
  r_post.resize( loci, n_subpops, last_pre + t_add -1);

  

  for (int i = last_pre-1; i < last_pre+t_add-1; i++){
    
    // Rcout << "\n\n\n";
    // Rcout << "i\n";
    // Rcout << i;
    // Rcout << "\n\n\n";
    // 
    // Rcout << "pAb_pre\n";
    // Rcout << pAb_pre;
    // Rcout << "\n";
    
    pAB_post.slice(i) = (arma::mat)pAB_pre.slice(i) * M;
    pAb_post.slice(i) = (arma::mat)pAb_pre.slice(i) * M;
    paB_post.slice(i) = (arma::mat)paB_pre.slice(i) * M;
    pab_post.slice(i) = (arma::mat)pab_pre.slice(i) * M;

    // Rcout << "M\n";
    // Rcout << M;
    // Rcout << "\n";
    // Rcout << "pAb_post\n";
    // Rcout << pAb_post;
    // Rcout << "\n";
    
    pA_post.slice(i) = pAB_post.slice(i) + pAb_post.slice(i);
    pB_post.slice(i) = pAB_post.slice(i) + paB_post.slice(i);
    D_post.slice(i) = (pAB_post.slice(i) % pab_post.slice(i)) - (pAb_post.slice(i) % paB_post.slice(i));
    r_post.slice(i) = ( D_post.slice(i) % D_post.slice(i) ) / (pA_post.slice(i) % (1-pA_post.slice(i)) % pB_post.slice(i)%(1-pB_post.slice(i)));


    arma::mat gameteAB = pAB_post.slice(i) - c*D_post.slice(i);
    arma::mat gameteAb = pAb_post.slice(i) + c*D_post.slice(i);
    arma::mat gameteaB = paB_post.slice(i) + c*D_post.slice(i);
    arma::mat gameteab = pab_post.slice(i) - c*D_post.slice(i);
    
    // Rcout << "D_post\n";
    // Rcout << D_post;
    // Rcout << '\n';
    // Rcout << "gameteAb\n";
    // Rcout << gameteAb;
    // Rcout << '\n';

    for (int subpop = 0; subpop < n_subpops; subpop++){
      int N_temp = N_v(subpop);
      arma::mat prob = arma::join_horiz(gameteAB.col(subpop), gameteAb.col(subpop), gameteaB.col(subpop), gameteab.col(subpop));
      
      arma::mat temp = gsl_mmm(2*N_temp, prob);
      
      
      // Rcout << "and allocate";
      // Rcout << temp;
      // Rcout << arma::size(pAB_pre);
      pAB_pre.subcube(0,subpop,i+1,   loci-1,subpop, i+1) = temp.col(0)/(2*N_temp);
      pAb_pre.subcube(0,subpop,i+1,   loci-1,subpop, i+1) = temp.col(1)/(2*N_temp);
      paB_pre.subcube(0,subpop,i+1,   loci-1,subpop, i+1) = temp.col(2)/(2*N_temp);
      pab_pre.subcube(0,subpop,i+1,   loci-1,subpop,i+1) = temp.col(3)/(2*N_temp);

//       pAb_pre[,i+1,subpop]<-temp[,2]/(2*N_temp)
//       paB_pre[,i+1,subpop]<-temp[,3]/(2*N_temp)
//       pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
//       pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
     }
    
//     # Calculate statistics

    pA_pre.slice(i+1) = pAB_pre.slice(i+1) + pAb_pre.slice(i+1);
    pB_pre.slice(i+1) = pAB_pre.slice(i+1) + paB_pre.slice(i+1);
    D_pre.slice(i+1) = ((pAB_pre.slice(i+1) % pab_pre.slice(i+1)) - (pAb_pre.slice(i+1) % paB_pre.slice(i+1)));
    r_pre.slice(i+1) = D_pre.slice(i+1)%D_pre.slice(i+1) / (pA_pre.slice(i+1)%(1-pA_pre.slice(i+1))%pB_pre.slice(i+1)%(1-pB_pre.slice(i+1)));

 
//     pA_pre[,i+1,] <- pAB_pre[,i+1,] + pAb_pre[,i+1,]
//     pB_pre[,i+1,] <- pAB_pre[,i+1,] + paB_pre[,i+1,]
//     D_pre[,i+1,] <- (pAB_pre[,i+1,] * pab_pre[,i+1,] - pAb_pre[,i+1,] * paB_pre[,i+1,])
//     r_pre[,i+1,] <- D_pre[,i+1,]**2 / (pA_pre[,i+1,]*(1-pA_pre[,i+1,])*pB_pre[,i+1,]*(1-pB_pre[,i+1,]))
     }
    int tcurr = last_pre+t_add;
    List o;
    o["pAB_pre"] = pAB_pre;
    o["pAb_pre"] = pAb_pre;
    o["paB_pre"] = paB_pre;
    o["pab_pre"] = pab_pre;
    o["pAB_post"] = pAB_post;
    o["pAb_post"] = pAb_post;
    o["paB_post"] = paB_post;
    o["pab_post"] = pab_post;
    o["D_pre"] = D_pre;
    o["D_post"] = D_post;
    o["pA_pre"] = pA_pre;
    o["pA_post"] = pA_post;
    o["pB_pre"] = pB_pre;
    o["pB_post"] = pB_post;
    o["r_pre"] = r_pre;
    o["r_post"] = r_post;
    o["loci"] = loci;
    o["t_tot"] = "nothing";
    o["n_subpops"] = n_subpops;
    o["t_curr"] = tcurr;
    return o;
}

// [[Rcpp::export]]
void test_slice(arma::cube C){
  for (int i = 0; i < 1000000; i++){
    C.slice(2) = C.slice(1)-0.000001;
  }
}

// [[Rcpp::export]]
void test_col(arma::cube C){
  for (int i = 0; i < 1000000; i++){
    C.col(2) = C.col(1)-0.000001;
  }
}


// [[Rcpp::export]]
Rcpp::List propegate_population_inner_C(arma::mat pAB_pre_first, arma::mat pAb_pre_first, arma::mat paB_pre_first, arma::mat pab_pre_first, NumericVector N_v, arma::mat M, double c, int loci, int n_subpops, int last_pre, int t_add, int age){
  arma::cube pAB_pre(loci, n_subpops, t_add+1, arma::fill::zeros);
  arma::cube pAb_pre(loci, n_subpops, t_add+1, arma::fill::zeros);
  arma::cube paB_pre(loci, n_subpops, t_add+1, arma::fill::zeros);
  arma::cube pab_pre(loci, n_subpops, t_add+1, arma::fill::zeros);

  arma::cube pAB_post(loci, n_subpops, t_add+1, arma::fill::zeros);
  arma::cube pAb_post(loci, n_subpops, t_add+1, arma::fill::zeros);
  arma::cube paB_post(loci, n_subpops, t_add+1, arma::fill::zeros);
  arma::cube pab_post(loci, n_subpops, t_add+1, arma::fill::zeros);
    
  arma::cube r_post(loci, n_subpops, t_add, arma::fill::zeros);
  arma::cube D_post(loci, n_subpops, t_add, arma::fill::zeros);
  arma::cube pA_post(loci, n_subpops, t_add, arma::fill::zeros);
  arma::cube pB_post(loci, n_subpops, t_add, arma::fill::zeros); 
    

  pAB_pre.slice(0) = pAB_pre_first;
  pAb_pre.slice(0) = pAb_pre_first; 
  paB_pre.slice(0) = paB_pre_first;
  pab_pre.slice(0) = pab_pre_first;
  
  
  for (int i = 0; i < t_add; i++){
    
    // Rcout << "\n\n\n";
    // Rcout << "i\n";
    // Rcout << i;
    // Rcout << "\n\n\n";
    // 
    // Rcout << "pAb_pre\n";
    // Rcout << pAb_pre;
    // Rcout << "\n";

    
    pAB_post.slice(i) = (arma::mat)pAB_pre.slice(i) * M;
    pAb_post.slice(i) = (arma::mat)pAb_pre.slice(i) * M;
    paB_post.slice(i) = (arma::mat)paB_pre.slice(i) * M;
    pab_post.slice(i) = (arma::mat)pab_pre.slice(i) * M;
    // Rcout << "M\n";
    // Rcout << M;
    // Rcout << "\n";
    // Rcout << "pAb_post\n";
    // Rcout << pAb_post;
    // Rcout << "\n";
    
    pA_post.slice(i) = pAB_post.slice(i) + pAb_post.slice(i);
    pB_post.slice(i) = pAB_post.slice(i) + paB_post.slice(i);
    D_post.slice(i) = (pAB_post.slice(i) % pab_post.slice(i)) - (pAb_post.slice(i) % paB_post.slice(i));
    r_post.slice(i) = ( D_post.slice(i) % D_post.slice(i) ) / (pA_post.slice(i) % (1-pA_post.slice(i)) % pB_post.slice(i)%(1-pB_post.slice(i)));

    arma::mat gameteAB = pAB_post.slice(i) - c*D_post.slice(i);
    arma::mat gameteAb = pAb_post.slice(i) + c*D_post.slice(i);
    arma::mat gameteaB = paB_post.slice(i) + c*D_post.slice(i);
    arma::mat gameteab = pab_post.slice(i) - c*D_post.slice(i);
    
    // Rcout << "D_post\n";
    // Rcout << D_post;
    // Rcout << '\n';
    // Rcout << "gameteAb\n";
    // Rcout << gameteAb;
    // Rcout << '\n';
    for (int subpop = 0; subpop < n_subpops; subpop++){
      int N_temp = N_v(subpop);
      
      arma::mat prob = arma::join_horiz(gameteAB.col(subpop), gameteAb.col(subpop), gameteaB.col(subpop), gameteab.col(subpop));
      
      arma::mat temp = gsl_mmm(2*N_temp, prob);
      
      // Rcout << "and allocate";
      // Rcout << temp;
      // Rcout << arma::size(pAB_pre);
      pAB_pre.subcube(0,subpop,i+1,   loci-1,subpop,i+1) = temp.col(0)/(2*N_temp);
      pAb_pre.subcube(0,subpop,i+1,   loci-1,subpop,i+1)  = temp.col(1)/(2*N_temp);
      paB_pre.subcube(0,subpop,i+1,   loci-1,subpop,i+1)  = temp.col(2)/(2*N_temp);
      pab_pre.subcube(0,subpop,i+1,   loci-1,subpop,i+1) = temp.col(3)/(2*N_temp);
      
      //       pAb_pre[,i+1,subpop]<-temp[,2]/(2*N_temp)
      //       paB_pre[,i+1,subpop]<-temp[,3]/(2*N_temp)
      //       pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
      //       pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
    }
    //     # Calculate statistics
    
    // pA_pre.col(i+1) = pAB_pre.col(i+1) + pAb_pre.col(i+1);
    // pB_pre.col(i+1) = pAB_pre.col(i+1) + paB_pre.col(i+1);
    // D_pre.col(i+1) = ((pAB_pre.col(i+1) % pab_pre.col(i+1)) - (pAb_pre.col(i+1) % paB_pre.col(i+1)));
    // r_pre.col(i+1) = D_pre.col(i+1)%D_pre.col(i+1) / (pA_pre.col(i+1)%(1-pA_pre.col(i+1))%pB_pre.col(i+1)%(1-pB_pre.col(i+1)));
    
    
    //     pA_pre[,i+1,] <- pAB_pre[,i+1,] + pAb_pre[,i+1,]
    //     pB_pre[,i+1,] <- pAB_pre[,i+1,] + paB_pre[,i+1,]
    //     D_pre[,i+1,] <- (pAB_pre[,i+1,] * pab_pre[,i+1,] - pAb_pre[,i+1,] * paB_pre[,i+1,])
    //     r_pre[,i+1,] <- D_pre[,i+1,]**2 / (pA_pre[,i+1,]*(1-pA_pre[,i+1,])*pB_pre[,i+1,]*(1-pB_pre[,i+1,]))
  }
  int t_curr = 1+t_add;
  List o;
  o["pAB_pre"] = pAB_pre;
  o["pAb_pre"] = pAb_pre;
  o["paB_pre"] = paB_pre;
  o["pab_pre"] = pab_pre;
  o["pAB_post"] = pAB_post;
  o["pAb_post"] = pAb_post;
  o["paB_post"] = paB_post;
  o["pab_post"] = pab_post;
  //o["D_pre"] = D_pre;
  o["D_post"] = D_post;
  //o["pA_pre"] = pA_pre;
  o["pA_post"] = pA_post;
  //o["pB_pre"] = pB_pre;
  o["pB_post"] = pB_post;
  //o["r_pre"] = r_pre;
  o["r_post"] = r_post;
  o["loci"] = loci;
  o["t_tot"] = "nothing";
  o["n_subpops"] = n_subpops;
  o["t_curr"] = t_curr;
  o["age"] = age+t_add;
  return o;
}




// [[Rcpp::export]]
Rcpp::List propegate_population_slim(List pop, arma::vec N_v, arma::mat M, double c, int t_add){
  
  
  int loci = as<int>(pop["loci"]);
  int t_min = as<int>(pop["t_curr"]);
  int n_subpops = as<int>(pop["n_subpops"]);
  
  arma::cube pAB_pre = as<arma::cube>(pop["pAB_pre"]);
  arma::cube pAb_pre = as<arma::cube>(pop["pAb_pre"]);
  arma::cube paB_pre = as<arma::cube>(pop["paB_pre"]);
  arma::cube pab_pre = as<arma::cube>(pop["pab_pre"]);
  
  arma::cube pAB_post = as<arma::cube>(pop["pAB_post"]);
  arma::cube pAb_post = as<arma::cube>(pop["pAb_post"]);
  arma::cube paB_post = as<arma::cube>(pop["paB_post"]);
  arma::cube pab_post = as<arma::cube>(pop["pab_post"]);
  
  arma::cube pA_pre = as<arma::cube>(pop["pA_pre"]);
  arma::cube pA_post = as<arma::cube>(pop["pA_post"]);
  
  arma::cube pB_pre = as<arma::cube>(pop["pB_pre"]);
  arma::cube pB_post = as<arma::cube>(pop["pB_post"]);
  
  arma::cube D_pre = as<arma::cube>(pop["D_pre"]);
  arma::cube D_post = as<arma::cube>(pop["D_post"]);
  arma::cube r_pre = as<arma::cube>(pop["r_pre"]);
  arma::cube r_post = as<arma::cube>(pop["r_post"]);
  
  
  
  //
  
  
  pAB_pre.resize( loci, t_min+1 + t_add, n_subpops );
  pAb_pre.resize( loci, t_min+1 + t_add, n_subpops );
  paB_pre.resize( loci, t_min+1 + t_add, n_subpops );
  pab_pre.resize( loci, t_min+1 + t_add, n_subpops );
  
  pAB_post.resize( loci, t_min+1 + t_add, n_subpops );
  pAb_post.resize( loci, t_min+1 + t_add, n_subpops );
  paB_post.resize( loci, t_min+1 + t_add, n_subpops );
  pab_post.resize( loci, t_min+1 + t_add, n_subpops );
  
  pA_pre.resize( loci, t_min+1 + t_add, n_subpops );
  pA_post.resize( loci, t_min+1 + t_add, n_subpops );
  
  pB_pre.resize( loci, t_min+1 + t_add, n_subpops );
  pB_post.resize( loci, t_min+1 + t_add, n_subpops );
  
  D_pre.resize( loci, t_min+1 + t_add, n_subpops );
  D_post.resize( loci, t_min+1 + t_add, n_subpops );
  r_pre.resize( loci, t_min+1 + t_add, n_subpops );
  r_post.resize( loci, t_min+1 + t_add, n_subpops );
  
  
  
  for (int i = t_min-1; i <= t_min+t_add-1; i++){
    

    
    pAB_post.col(i) = (arma::mat)pAB_pre.col(i) * M;
    pAb_post.col(i) = (arma::mat)pAb_pre.col(i) * M;
    paB_post.col(i) = (arma::mat)paB_pre.col(i) * M;
    pab_post.col(i) = (arma::mat)pab_pre.col(i) * M;
    

    
    pA_post.col(i) = pAB_post.col(i) + pAb_post.col(i);
    pB_post.col(i) = pAB_post.col(i) + paB_post.col(i);
    D_post.col(i) = (pAB_post.col(i) % pab_post.col(i) - pAb_post.col(i) % paB_post.col(i));
    r_post.col(i) = ( D_post.col(i) % D_post.col(i) ) / (pA_post.col(i) % (1-pA_post.col(i)) % pB_post.col(i)%(1-pB_post.col(i)));
    
    
    arma::mat gameteAB = pAB_post.col(i) - c*D_post.col(i);
    arma::mat gameteAb = pAb_post.col(i) + c*D_post.col(i);
    arma::mat gameteaB = paB_post.col(i) + c*D_post.col(i);
    arma::mat gameteab = pab_post.col(i) + c*D_post.col(i);
    

    
    for (int subpop = 0; subpop < n_subpops; subpop++){
      int N_temp = N_v(subpop);
      arma::mat prob = arma::join_horiz(gameteAB.col(subpop), gameteAb.col(subpop), gameteaB.col(subpop), gameteab.col(subpop));
      //Rcout << prob;
      arma::mat temp = gsl_mmm(2*N_temp, prob);
      // Rcout << "and allocate";
      // Rcout << temp;
      // Rcout << arma::size(pAB_pre);
      pAB_pre.subcube(0,i+1,subpop,loci-1,i+1,subpop) = temp.col(0)/(2*N_temp);
      pAb_pre.subcube(0,i+1,subpop,loci-1,i+1,subpop) = temp.col(1)/(2*N_temp);
      paB_pre.subcube(0,i+1,subpop,loci-1,i+1,subpop) = temp.col(2)/(2*N_temp);
      pab_pre.subcube(0,i+1,subpop,loci-1,i+1,subpop) = temp.col(3)/(2*N_temp);
      
      //       pAb_pre[,i+1,subpop]<-temp[,2]/(2*N_temp)
      //       paB_pre[,i+1,subpop]<-temp[,3]/(2*N_temp)
      //       pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
      //       pab_pre[,i+1,subpop]<-temp[,4]/(2*N_temp)
    }
    
    //     # Calculate statistics
    
    
    
    
    //     pA_pre[,i+1,] <- pAB_pre[,i+1,] + pAb_pre[,i+1,]
    //     pB_pre[,i+1,] <- pAB_pre[,i+1,] + paB_pre[,i+1,]
    //     D_pre[,i+1,] <- (pAB_pre[,i+1,] * pab_pre[,i+1,] - pAb_pre[,i+1,] * paB_pre[,i+1,])
    //     r_pre[,i+1,] <- D_pre[,i+1,]**2 / (pA_pre[,i+1,]*(1-pA_pre[,i+1,])*pB_pre[,i+1,]*(1-pB_pre[,i+1,]))
  }
  int tcurr = t_min+t_add;
  List o;
  o["pAB_pre"] = pAB_pre;
  o["pAb_pre"] = pAb_pre;
  o["paB_pre"] = paB_pre;
  o["pab_pre"] = pab_pre;
  o["pAB_post"] = pAB_post;
  o["pAb_post"] = pAb_post;
  o["paB_post"] = paB_post;
  o["pab_post"] = pab_post;
  o["D_post"] = D_post;
  o["pA_post"] = pA_post;
  o["pB_post"] = pB_post;
  o["loci"] = loci;
  o["t_tot"] = "nothing";
  o["n_subpops"] = n_subpops;
  o["t_curr"] = tcurr;
  return o;
  
  
}





