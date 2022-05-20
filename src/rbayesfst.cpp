#include <algorithm>
#include <ctime>
#include <random>
#include <Rcpp.h>
//#include <gperftools/profiler.h> // only needed for profiling

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;
using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


class BayesFst{
    /* USER-DEFINED CONSTANTS */
    /* PRIORS FOR a_i, b_j AND g_ij - NORMAL DISTRIBUTION (mu , si^2)*/
    
private:
    double alphaMu; // (expected) mean of Normal distribution for prior of a_i
    double alphaSigma; // (expected) sd of Normal distribution for prior of a_i
    double betaMu; // mean of Normal distribution for prior of b_j
    double betaSigma; // sd of Normal distribution for prior of b_j
    double gammaMu; // mean of Normal distribution for prior of g_ij
    double gammaSigma; // sd of Normal distribution for prior of g_ij
    
    /* CONSTANT FOR PROPOSAL DISTRIBUTION */
    
    double uSigma; // scale param for normal updates */
    double pSigma; // scale param for Dirichlet updates of p_ij */
    
    /* CORRELATION PARAMETER (FIXED) BETWEEN ADJACENT MARKERS */
    
    double cor; // the value of cor determines the shrinkage value
    // associated with the sd of each a_i so that the
    // overall - expected - sd is alphaSigma (defined above).
    // The shrinkage value (srk) is defined below
    
    double srk; // shrinkage value for sd of each a_i (see above)
    
    // input member variables
    
    int m_nLoci, m_nPops;
    unsigned int m_nSeed;
    long m_nKeep, m_nAcceptanceRateGap, m_nNumIt, m_nDiscard;
    double rtp, ftp;
    bool m_bGammaSwitch;
    IntegerVector m_numAlleles;
    vector< vector <int> > m_popSums;
    vector< vector< vector<int> > > m_Counts;
    bool m_bDataLoaded;
    
    // input variables which were constants and so may be rarely changed
    int m_nNumOut; // default 2001
    
    // Random number generator stuff
    mt19937 m_mtEngine;
    normal_distribution<double> m_rngZ;
    uniform_real_distribution<double> m_rngU;
    
    // for the Wichman-Hill RNG - need for debugging but will go eventually
    long ix, iy, iz;
    
    // internal variables
    
    vector<double> alpha, beta, gamma;
    vector<double> lp0, lp1;
    vector< vector<double> > pi;
    bool bIllegal;
    vector<double> m_vdFst;
    vector<int> m_vnNumZeros;
    double m_dMeanFst;
    double usa, usb, usg;

    double calcFst(double w){
        /* calculate Fst(loc_i , pop_j) from hyperparameters a_i, b_j and g_ij */
        
        return std::exp(w) / (1 + std::exp(w));
    }
    
    void initFst(){
        // generate the initial Fst values and sort them
        double mu = alphaMu + betaMu + (m_bGammaSwitch ? gammaMu : 0);
        double sigma = std::sqrt(alphaSigma * alphaSigma + betaSigma * betaSigma + (m_bGammaSwitch ? gammaSigma * gammaSigma: 0));
        m_dMeanFst = 0;
        
        for(int i = 0; i < m_nNumOut; i++){
            m_vdFst[i] =  calcFst(rnorm(mu, sigma));
            m_dMeanFst += m_vdFst[i];
        }
        
        std::sort(m_vdFst.begin(), m_vdFst.end());
        m_dMeanFst /= m_nNumOut;
    }
    
    void initAlleleFreqs(){ // initialize the allele frequencies
        pi.resize(m_nLoci);
        
        for(int loc = 0; loc < m_nLoci; loc++){
            
            pi[loc].resize(m_numAlleles[loc]);
            fill(pi[loc].begin(), pi[loc].end(), 0);
            
            for(int a = 0; a < m_numAlleles[loc]; a++){
                for(int pop = 0; pop < m_nPops; pop++){
                    if(m_popSums[loc][pop] > 0){
                        pi[loc][a] += (double)(m_Counts[loc][pop][a] + 1) / (double)((m_popSums[loc][pop] + m_numAlleles[loc]) * (m_nPops - m_vnNumZeros[loc]));
                    }
                }
            }
        }
    }
    
    void init_gen(unsigned int seed){
        m_mtEngine = mt19937(seed);
        m_rngU = uniform_real_distribution<>(0.0, 1.0);
        m_rngZ = normal_distribution<double>(0.0, 1.0);
    }
    
    void initWichHill(long x, long y, long z){
      ix = x;
      iy = y;
      iz = z;
    }
    
    void initHyperparams(){
        alpha.resize(m_nLoci);
        
        alpha[0] = rnorm(alphaMu, alphaSigma);
        for(int i = 1; i < m_nLoci; i++)
            alpha[i] = rnorm(cor * alpha[i-1], srk * alphaSigma);
        
        beta.resize(m_nPops);
        
        if(m_bGammaSwitch)
            gamma.resize(m_nPops * m_nLoci);
        
        for(int j = 0; j < m_nPops; j++){
            beta[j] = rnorm(betaMu, betaSigma);
            
            if(m_bGammaSwitch){
                for(int i = 0; i < m_nLoci; i++){
                    gamma[i * m_nPops + j] = rnorm(gammaMu, gammaSigma);
                }
            }
        }
    }
    
    double initLogPostDens(const vector<double>& a, const vector<double>& b, const vector<double>& g ,
                           vector<double>& lp, vector<double>& pa, vector<double>& pg){
        double lpd = 0;
        
        for(int j = 0; j < m_nPops; j++)
            lpd += logdnorm(b[j] , betaMu , betaSigma);     /* PRIOR FOR b_j */
        
        for(int i = 0;i < m_nLoci; i++){
            if(i==0){
                lp[i] = logdnorm(a[i], alphaMu, alphaSigma);	/* PRIOR FOR a_i */
            }else if(i > 0 && i < m_nLoci){
                lp[i] = logdnorm(a[i], cor * a[i-1], srk * alphaSigma);
                pa[i] = lp[i];
            }

            pa[m_nLoci] = 0;
            
            if(m_bGammaSwitch){
              for (int j = 0; j < m_nPops; j++){
                pg[i] += logdnorm(g[i], gammaMu, gammaSigma);
              }
              lp[i] += pg[i];
            }

            for (int j = 0; j < m_nPops; j++){
                if (m_popSums[i][j] > 0){
                    lp[i] += (m_bGammaSwitch ? logMultinomDirichlet( calcFst(a[i] + b[j] + g[i * m_nPops + j]), i, j) : logMultinomDirichlet( calcFst(a[i] + b[j]), i, j));   /* LIKELIHOOD */
                }
            }
            lpd += lp[i];
        }
        
        return lpd;
    }
    
    double logdnorm(double lx, double mu, double si){
        /* evaluate (Normal) prior log densities for the a_i, b_j and g_ij
         (excluding constant) */
        
        return -0.5 * std::pow((lx - mu) / si, 2);
    }
    
    double logDirichlet(const vector<double>& par, const vector<double>& vec){
        double sum;
        double konst;
        
        sum = 0.0;
        konst = std::accumulate(par.begin(), par.end(), 0.0);
        
        for(int j = 0; j < (int)par.size(); j++){
            sum += (par[j] - 1.0) * std::log(vec[j]);
            //if(debug) Rprintf("j = %d, sum = %f, par[j] = %f, vec[j] = %f\n", j, sum, par[j], vec[j]);
            sum -= std::lgamma(par[j]);
            //if(debug) Rprintf("%f\n", sum);
        }
        //if(debug) Rprintf("%f %f %f\n", sum, konst, sum + lgamma(konst));
        return sum + std::lgamma(konst);
    }
    
    
    double logMultinomDirichlet(double f0, int i, int j){
        /* modified version of function by M.A. Beaumont for evaluating
         multinomial-Dirichlet likelihood for locus i, population j
         (excluding constant)
         */
        
        double t1 = 0.0, a;
        double r = (1 - f0) / f0;
        
        for (int k = 0; k< m_numAlleles[i]; k++){
            a = pi[i][k] * r;
            t1 += std::lgamma(m_Counts[i][j][k] + a) -std::lgamma(a);
        }
        t1 += std::lgamma(r) - std::lgamma(m_popSums[i][j] + r);
        return t1;
    }
    
    vector<double> pargen(const vector<double>& oldpi){
        vector<double> fpar, rpar, tvec;
        int k;
        double cum = 0.0;
        int nA = oldpi.size();
        vector<double> newpi;
        
        newpi.reserve(nA);
        fpar.reserve(nA);
        rpar.reserve(nA);
        tvec.reserve(nA);
        
        for(k = 0; k < nA; k++){
            tvec.push_back(oldpi[k]);
            
            if(tvec[k] < 1.0e-3){
                tvec[k] = 1.0e-3;
            }
            
            fpar.push_back(pSigma * tvec[k]);
            newpi.push_back(::Rf_rgamma(fpar[k], 1.0)); // rgamma is from Rcpp so it is vectorized
            
            if(newpi[k] == 0.0){
                bIllegal = 1;
                return newpi;
            }
            
            cum += newpi[k];
        }
        
        for(k = 0; k < nA; k++){
            newpi[k] /= cum;
            tvec[k] = newpi[k];
            
            if(tvec[k] < 1.0e-3){
                tvec[k] = 1.0e-3;
            }
            rpar.push_back(pSigma * tvec[k]);
        }
        
        /* forward and reverse  log-Hastings terms */
        ftp = logDirichlet(fpar, newpi);
        rtp = logDirichlet(rpar, oldpi);
        
        return newpi;
    }
    
    void printHyperparams(void){
      for(int i = 0; i < m_nLoci; i++)
        Rprintf("alpha[%d] = %6.2f\n", i + 1, alpha[i]);
      
      for(int j = 0; j < m_nPops; j++)
        Rprintf("beta[%d] = %6.2f\n", j + 1, beta[j]);
      
      if(m_bGammaSwitch){
        for(int i = 0; i < m_nLoci; i++){
          for(int j = 0; j < m_nPops; j++){
            Rprintf("%6.2f",gamma[i * m_nPops + j]);
          }
          Rprintf("\n");
        }
      }
    }
    
    
    double quantileFst(double p){
        double idx = p * m_nNumOut + 0.5;
        if(idx - std::floor(idx) < 0.1){
            return m_vdFst[(int)idx];
        } //else{
        int k = (int)idx;
        return m_vdFst[k] + fmod(idx, 1.0) * (m_vdFst[k + 1] - m_vdFst[k]);
    }
    
    double WichHill(){
      // Wichmann & Hill random number generator (Algorithm AS183 in Applied
      // Statistics 31 (1982)) 
      //
      
      ix = (171 * ix) % 30269;
      iy = (172 * iy) % 30307;
      iz = (170 * iz) % 30323;
      return fmod((ix / 30269.0) + (iy / 30307.0) + (iz / 30323.0), 1.0);
    }
    
    
    double runif(void){
        return m_rngU(m_mtEngine);
    }
    
    double rnorm(double mu = 0, double sigma = 1){
        return sigma * m_rngZ(m_mtEngine) + mu;
    }
    
    // for debugging - eventually will be removed
    double Rnorm(double mu = 0, double sigma = 1){
      /* Algorithm 3.6 (Polar) in Ripley (1987) for simulating normal RVs */
      
      double u1, u2, w;
      
      do{
        u1 = 2 * WichHill() - 1;
        u2 = 2 * WichHill() - 1;
        w = u1 * u1 + u2 * u2;
      }while(w <= 0 || w >= 1);
      
      double z = u1 * std::sqrt(-2 * std::log(w) / w);
      return sigma * z + mu;
    }
    
    int update_alpha(vector<double>&a, const vector<double>& b, vector<double>& g,
                 double& lpd, vector<double>& lp,  vector<double>& pa, vector<double>& pg){
        
        int i, j;
        vector<double> alphaDash, lp1, pa1;
        vector<double> gammaDash, pg1;
        int nAccept = 0;
        
        alphaDash.resize(m_nLoci);
        lp1.resize(m_nLoci);
        pa1.resize(m_nLoci);
        
        if(m_bGammaSwitch){
            pg1.resize(m_nLoci);
            gammaDash.resize(m_nPops * m_nLoci);
            fill(pg1.begin(), pg1.end(), 0);
        }
        
        for (i = 0; i < m_nLoci; i++){
            vector<double>  oldpi(pi[i]);
            
            // choose new p and new alpha
            pi[i] = pargen(oldpi);
            alphaDash[i] = rnorm(a[i], usa);
            
            if(i == 0){
                lp1[0] = logdnorm(alphaDash[0] , alphaMu , alphaSigma);
                pa1[0] = logdnorm(alphaDash[0] , alphaMu , alphaSigma);
                pa1[1] = logdnorm(a[1] , cor * alphaDash[0] , srk * alphaSigma);
            }
            else if(i > 0 &&i < m_nLoci - 1){
                lp1[i] = logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i] = logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i + 1]= logdnorm(a[i + 1] , cor * alphaDash[i] , srk * alphaSigma);
            }
            else if(i == m_nLoci - 1){
                lp1[i] = logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i] = logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i + 1] = 0;
            }
            
            if(m_bGammaSwitch){
                for (int j = 0; j < m_nPops; j++){
                    gammaDash[i * m_nPops + j] = rnorm(g[i * m_nPops + j], usg);
                    pg1[i] += logdnorm(gammaDash[i * m_nPops + j], gammaMu, gammaSigma);
                }
            }
            
            /* EVALUATE LOG - POSTERIOR DENSITY FOR CANDIDATE VECTORS (newz) */
            
            for (j = 0; j < m_nPops; j++){
                if (m_popSums[i][j] > 0){
                    lp1[i] += (m_bGammaSwitch ? logMultinomDirichlet(calcFst(alphaDash[i] + b[j] + gammaDash[i * m_nPops + j]) , i, j) : logMultinomDirichlet(calcFst(alphaDash[i] + b[j]) , i, j));
                }
            }
            
            /* DECIDE WHETHER TO ACCEPT/REJECT CANDIDATE VECTORS */
            
            double logU = std::log(runif());
           // Rprintf("%f %f %f %f %f %f %f\n", lpd, lp1[i], pa1[i + 1], lp[i], pa[i + 1], rtp, ftp);
            
            if(!bIllegal && (logU < lp1[i] + pa1[i + 1] - lp[i] - pa[i + 1] + rtp - ftp)){/* ACCEPT */
                lpd += lp1[i] + pa1[i + 1] - lp[i] - pa[i + 1];
                lp[i] = lp1[i];
                lp[i + 1] += pa1[i + 1] - pa[i + 1];
                pa[i] = pa1[i];
                pa[i + 1] = pa1[i + 1];
                a[i] = alphaDash[i];
                if(m_bGammaSwitch){
                    pg[i] = pg1[i];
                    for (j = 0; j < m_nPops; j++){
                        g[i * m_nPops + j] = gammaDash[i * m_nPops + j];
                    }
                }
                nAccept++;
            }
            else{
                pi[i] = oldpi;
                
                if(i > 0)
                    pa[i]= logdnorm(a[i] , cor * a[i - 1] , srk * alphaSigma);
            }
        }
        
        return nAccept;
    }
    
    int update_beta(const vector<double>&a, vector<double>& b, const vector<double>& g,
                 double& lpd, vector<double>& lp, const vector<double>& pa){
        
        /* CHOOSE betaDash FOR POTENTIAL UPDATE */
        vector<double> betaDash;
        betaDash.resize(b.size());
        
        
        
        for (int j = 0; j < m_nPops; j++)
            betaDash[j] = rnorm(b[j], usb);
        
        /* EVALUATE LOG-POSTERIOR DENSITY FOR CANDIDATE VECTORS (newz) */
        
        double lpd1 = 0;
        
        for(int j = 0; j < m_nPops; j++)
            lpd1 += logdnorm(betaDash[j], betaMu, betaSigma);
        
        vector<double> lp1;
        lp1.resize(m_nLoci);
        
        for (int i = 0; i < m_nLoci; i++){
            lp1[i] = pa[i];

            for (int j = 0; j < m_nPops; j++){

                if (m_popSums[i][j] > 0){
                    lp1[i]+= (m_bGammaSwitch ? logMultinomDirichlet(calcFst(a[i] + betaDash[j] + g[i * m_nPops + j]), i, j) : logMultinomDirichlet(calcFst(a[i] + betaDash[j]), i, j));
                }
            }
            lpd1 += lp1[i];
        }
        
       // Rprintf("%f %f\n", lpd, lpd1);

        /* DECIDE WHETHER TO ACCEPT/REJECT CANDIDATE VECTORS */
        double logU = std::log(runif());

        if (!bIllegal && logU < lpd1 - lpd){ // ACCEPT
            lpd = lpd1;
            lp = lp1;
            b = betaDash;
            return 1;
        }//else
        
        return 0;
    }
    
public:
    BayesFst(void){
        // default constructor
        alphaMu = 0.0;     /* (expected) mean of Normal distribution for prior of a_i */
        alphaSigma = 1.0;     /* (expected) sd of Normal distribution for prior of a_i */
        betaMu = -2.0;     /* mean of Normal distribution for prior of b_j */
        betaSigma = 1.8;     /* sd of Normal distribution for prior of b_j */
        gammaMu = 0.0;     /*  mean of Normal distribution for prior of g_ij */
        gammaSigma = 0.5;     /* sd of Normal distribution for prior of g_ij */
        
        /* CONSTANT FOR PROPOSAL DISTRIBUTION */
        
        uSigma = 0.5;     /* scale param for normal updates */
        pSigma = 1000;     /* scale param for Dirichlet updates of p_ij */
        
        //
        usa = alphaSigma * uSigma;
        usb = betaSigma * uSigma / std::sqrt(m_nLoci);
        usg = gammaSigma * uSigma;
        
        /* CORRELATION PARAMETER (FIXED) BETWEEN ADJACENT MARKERS */
        
        cor = 0; // the value of cor determines the shrinkage value
        // associated with the sd of each a_i so that the
        // overall - expected - sd is alphaSigma (defined above).
        // The shrinkage value (srk) is defined below
        
        srk = 1;
        
        // Interaction is switched OFF by default
        m_bGammaSwitch = false;
        
        // Data is not loaded by default
        m_bDataLoaded = false;
        
        // Set the number of iteration to -1, so that the system can detect
        // whether the user has changed this (to something sensible)
        m_nNumIt = -1;
        
    }
  
    
    
    bool getInteraction(){
        return m_bGammaSwitch;
    }
    
    void setInteraction(bool bGammaSwitch){
        m_bGammaSwitch = bGammaSwitch;
    }
    
    void printData(void){
        Rprintf("This data set has:\n");
        Rprintf("%d Populations:\n", m_nPops);
        Rprintf("%d Loci:\n", m_nLoci);
        for(int i = 0; i < m_nLoci; i++)
            Rprintf("Locus %d has %d alleles\n", i + 1, m_numAlleles[i]);
    }
  
    void printCounts(void){
      Rprintf("Loc  Pop  Sample size  Allele counts\n");
      
      for(int loc = 0; loc < m_nLoci; loc++){
        for(int pop = 0; pop < m_nPops; pop++){
          Rprintf("%4d%5d%5d", loc + 1, pop + 1, m_popSums[loc][pop]);
          for(int a  = 0; a < m_numAlleles[loc]; a++){
            Rprintf("%8d", m_Counts[loc][pop][a]);
          }
          Rprintf("\n");
        }
      }
    }
  
  void printInitialPvals(void){
    Rprintf("****** Initial p values ******\n\n");
    
    for (int i = 0; i < m_nLoci; i++){
      Rprintf("Locus %4d", i + 1);
      for (int k = 0; k < m_numAlleles[i]; k++){ /* LAPLACE VALUE FOR ALLELE FREQS */
          Rprintf("%7.3f", pi[i][k]);
      }
      Rprintf("\n");
    }
    
    Rprintf("\n");
    
  }
  
  
    void printRunInfo(void){
      Rprintf("****** Settings for MCMC algorithm ******\n\n");
      
      if(m_nNumIt > 0){
        Rprintf("Burn-in length (iterations): %ld\n", m_nDiscard);
        Rprintf("Thinning interval: %ld\n", m_nKeep);
        Rprintf("Number of iterations (post burn-in): %ld\n" , m_nNumIt);
        Rprintf("Number of outputs: %ld\n", m_nNumOut);
      }else{
        Rprintf("Run parameters have not been set yet\n");
      }
    }
    
    void printFstSummary(){
        Rprintf("\n** Prior summary stats, based on %ld simulated values, for each F_ST(i,j) **\n\n",m_nNumOut);
        
        Rprintf("-----------------------------------------------------------------------------------\n");
        Rprintf("|F_ST(i,j)|Mean    |2.5%c-ile| 5 %c-ile|25 %c-ile| median |75 %c-ile|95 %c-ile|97.5%cile|\n",37,37,37,37,37,37);
        Rprintf("|---------------------------------------------------------------------------------|\n");
        Rprintf("|Prior    |%7f|", m_dMeanFst);
        
        vector<double> vdPercentiles = {0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975};
        
        for(std::vector<double>::iterator it = vdPercentiles.begin(); it != vdPercentiles.end(); ++it){
            Rprintf("%8f|", quantileFst(*it));
        }
        
        Rprintf("\n");
        
    }
  
    List getPriorParameters(void){
      List l;
      
      l["alpha"] = NumericVector::create(alphaMu, alphaSigma);
      l["beta"] = NumericVector::create(betaMu, betaSigma);
      l["gamma"] = NumericVector::create(gammaMu, gammaSigma);
      l["uSigma"] = uSigma;
      l["pSigma"] = pSigma;
      l["cor"] = cor;
      
      return l;
      
    }
    
    void setPriorParameters(NumericVector alphaPrior, NumericVector betaPrior, NumericVector gammaPrior,
                            double uSigma, double pSigma, double shrink){
        alphaMu = alphaPrior[0];
        alphaSigma = alphaPrior[1];
        betaMu = betaPrior[0];
        betaSigma = betaPrior[1];
        gammaMu = gammaPrior[0];
        gammaSigma = gammaPrior[1];
        
        this->uSigma = uSigma;
        this->pSigma = pSigma;
        
        cor = shrink;
        srk = std::sqrt(1 - cor * cor);
        
        usa = alphaSigma * uSigma;
        usb = betaSigma * uSigma / std::sqrt(m_nLoci);
        usg = gammaSigma * uSigma;
    }
    
    void setData(List data){
        m_nPops = data["nPops"];
        m_nLoci = data["nLoci"];
        
        usb = betaSigma * uSigma / std::sqrt(m_nLoci);
        
        m_numAlleles = as<vector<int> >(data["numAlleles"]);
        
        NumericMatrix popSums = as<NumericMatrix>(data["locusPopSums"]);
        m_popSums.resize(m_nLoci);
        
        for(int i = 0; i < m_nLoci; i++){
            NumericVector r = popSums(i, _);
            m_popSums[i] = as< vector<int> >(r);
        }
        
        List listCounts = as<List>(data["dbCounts"]);
        m_Counts.resize(m_nLoci);
        m_vnNumZeros.resize(m_nLoci);
        
        fill(m_vnNumZeros.begin(), m_vnNumZeros.end(), 0);
        
        for(int loc = 0; loc < m_nLoci; loc++){
            
            m_Counts[loc].resize(m_nPops);
            IntegerMatrix locusCounts = as<IntegerMatrix>(listCounts[loc]);
            
            for(int pop = 0; pop < m_nPops; pop++){
                
                if(m_popSums[loc][pop] == 0){
                    m_vnNumZeros[loc]++;
                }
                
                m_Counts[loc][pop].resize(m_numAlleles[loc]);
                
                
                for(int a = 0; a < m_numAlleles[loc]; a++){
                    m_Counts[loc][pop][a] = locusCounts(pop, a);
                }
            }
        }
        
        m_bDataLoaded = true;
        Rprintf("Input data successfully initialized\n");
    }
  
    bool isDataLoaded(void){
      return m_bDataLoaded;
    }
    
    void setRunParameters(int numOut = 21, double keepPpn = 0.2, double discardPpn = 0.05, double acceptPpn = 0.02, bool bPrint = false){
        m_nNumOut = numOut;
        m_nKeep = (long)(keepPpn * m_nPops * m_nLoci);      /* interval between successive outputs */
        m_nNumIt = m_nNumOut * m_nKeep; /* number of iterations of alg. after burn-n */
        m_nDiscard = (long)(discardPpn * m_nNumIt);   /* length of burn-in  */
        m_nAcceptanceRateGap = (long)(acceptPpn * (m_nNumIt + m_nDiscard));  /* gap between output of acceptance rates */
        
        m_vdFst.resize(m_nNumOut);
        
        if(bPrint){
            Rprintf("****** Settings for MCMC algorithm ******\n\n");
            Rprintf("Burn-in length (iterations): %ld\n", m_nDiscard);
            Rprintf("Thinning interval: %ld\n", m_nKeep);
            Rprintf("Number of iterations (post burn-in): %ld\n", m_nNumIt);
            Rprintf("Number of outputs: %ld\n", m_nNumOut);
            
            Rprintf("Using Normal prior for hyperparameters alpha_i with mean %.2f, sd %.2f\n",alphaMu,alphaSigma);
            Rprintf("Using Normal prior for hyperparameters beta_j with mean %.2f, sd %.2f\n",betaMu,betaSigma);
            
            if (m_bGammaSwitch){
                Rprintf("Using interaction terms in model\n");
                Rprintf("Using Normal prior for hyperparameters gamma_ij with mean %.2f, sd %.2f\n\n",gammaMu,gammaSigma);
            }
        }
    }
    
    List run(unsigned int seed){
        m_nSeed = seed;
        init_gen(m_nSeed);

        initFst();
        printFstSummary();
        
        /* INITIALISE HYPERPARAMETER VECTORS (a, b) - GENERATE FROM PRIOR */
        initHyperparams();
        printHyperparams();
        
        // Initialise allele frequencies
        initAlleleFreqs();
        printInitialPvals();
        
        /* INITIALISE LOG-POSTERIOR DENSITY VALUE (z) */

        vector<double> lp, pa, pg;
        lp.resize(m_nLoci);
        pa.resize(m_nLoci + 1);

        if(m_bGammaSwitch){
            pg.resize(m_nLoci);
            fill(pg.begin(), pg.end(), 0);
        }

        double logPosteriorDensity = initLogPostDens(alpha, beta, gamma, lp, pa, pg);
        //Rprintf("%f\n", logPosteriorDensity);


        /* RUN METROPOLIS ALGORITHM (for numit+discard iterations) */

        int jmp1 = 0;
        int jmp2 = 0;

        NumericMatrix postAlpha(m_nNumOut, m_nLoci);
        NumericMatrix postBeta(m_nNumOut, m_nPops);
        NumericMatrix postGamma(m_nNumOut, m_nLoci * m_nPops);
        NumericVector lpd(m_nNumOut);

        int sumAlleles = std::accumulate(m_numAlleles.begin(), m_numAlleles.end(), 0.0);
        NumericMatrix postP(m_nNumOut, sumAlleles);
        
        //clock_t begin = clock();
        //ProfilerStart("/Users/jcur002/Dropbox/Code/git/rbayesfst.log"); Only needed for profiling
        
        int ctr = 0;
        Progress p1(m_nDiscard, true);
        Progress p2(m_nNumIt, true);

        for (int nCurrentIteration = -m_nDiscard; nCurrentIteration <= m_nNumIt; nCurrentIteration++){
            bIllegal = false;

            jmp1 += update_beta(alpha, beta, gamma, logPosteriorDensity, lp, pa);
            jmp2 += update_alpha(alpha, beta, gamma, logPosteriorDensity, lp, pa, pg);

            if(nCurrentIteration < 0)
              p1.increment();
            else
              p2.increment();

            if(nCurrentIteration > 0){
                if ((nCurrentIteration / m_nKeep) * m_nKeep == nCurrentIteration){
                    lpd[ctr] = logPosteriorDensity;
                    postAlpha(ctr, _) = NumericVector(alpha.begin(), alpha.end());
                    postBeta(ctr, _) = NumericVector(beta.begin(), beta.end());

                    if(m_bGammaSwitch){
                      postGamma(ctr, _) = NumericVector(gamma.begin(), gamma.end());
                    }

                    int pos = 0;
                    for(int loc = 0; loc < m_nLoci; loc++){
                      for(int a = 0; a < m_numAlleles[loc]; a++){
                        postP(ctr, pos++) = pi[loc][a];
                      }
                    }

                    ctr++;
                }
            }
            if (((nCurrentIteration + m_nDiscard) > 0) && ((nCurrentIteration+m_nDiscard)/m_nAcceptanceRateGap)*m_nAcceptanceRateGap==(nCurrentIteration+m_nDiscard)){
              double percentDone = 100.0 * nCurrentIteration / (double) m_nNumIt;
              
              // see if the user wants to halt computation
              // Rcpp::checkUserInterrupt();
              if(nCurrentIteration < 0){
                if(p1.check_abort()){
                  Rprintf("User interrupted\n");
                  List emptyList;
                  return emptyList;
                }
              }else{
                if(p2.check_abort()){
                  Rprintf("User interrupted\n");
                  List emptyList;
                  return emptyList;
                }
              }
              
              //if(percentDone > 0)
              //  Rprintf("Percent: %5.2f Iter: %8d, Accpt. rate: %12.6f%12.6f\n", percentDone, nCurrentIteration, ((double)(jmp1)/m_nAcceptanceRateGap), ((double)(jmp2)/m_nAcceptanceRateGap/m_nLoci));
              
              jmp1=0;
              jmp2=0;
            }
        }
        /*ProfilerStop();
        
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        Rprintf("%f seconds elapsed\n", elapsed_secs);*/  // only needed for profiling

        List results;

        results["lpd"] = lpd;
        results["alpha"] = postAlpha;
        results["beta"] = postBeta;
        results["p"] = postP;
        if(m_bGammaSwitch)
          results["gamma"] = postGamma;
        results["nout"] = m_nNumOut;
        results["nloci"] = m_nLoci;
        results["npop"] = m_nPops;
        results["interaction"] = m_bGammaSwitch;

        return results;
    }
  
  double ldiriTest(NumericVector par, NumericVector vec){
    return logDirichlet(as< vector<double> >(par), as< vector<double> >(vec));
  }
  
  int getNumIterations(void){
    return m_nNumIt;
  }
  
  long getThin(void){
    return m_nKeep;
  }
};
    
    
    
RCPP_MODULE(BayesFst) {
    using namespace Rcpp;
    
    class_<BayesFst>( "BayesFst")
    .default_constructor("Standard constructor")
    .method("printData", &BayesFst::printData)
    .method("printCounts", &BayesFst::printCounts)
    .method("printInitialPvals", &BayesFst::printInitialPvals)
    .method("printFstSummary", &BayesFst::printFstSummary)
    .method("printRunInfo", &BayesFst::printRunInfo)
    .method("run", &BayesFst::run)
    .method("setData", &BayesFst::setData)
    .method("setRunParameters", &BayesFst::setRunParameters)
    .method("ldiriTest", &BayesFst::ldiriTest)
    .method("isDataLoaded", &BayesFst::isDataLoaded)
    .method("getPriorParams", &BayesFst::getPriorParameters)
    .method("setPriorParams", &BayesFst::setPriorParameters)
    .method("getNumIterations", &BayesFst::getNumIterations)
    .method("getThin", &BayesFst::getThin)
    .property("interaction", &BayesFst::getInteraction, &BayesFst::setInteraction)
    ;
}

//' @title Calculate Fst using the model 
//' @param results - the output of \code{\link{sample.bayesFst}}
//' @return a matrix with \code{nout} rows and \code{nloci * npop} columns
//' @export
// [[Rcpp::export]]
NumericMatrix CalcFst(List results){
  int nout = as<int>(results["nout"]);
  int nloc = as<int>(results["nloci"]);
  int npop = as<int>(results["npop"]);
  bool bInteraction = as<bool>(results["interaction"]);
  
  NumericMatrix fst(nout, nloc * npop);
  NumericMatrix alpha = as<NumericMatrix>(results["alpha"]);
  NumericMatrix beta = as<NumericMatrix>(results["beta"]);
  NumericMatrix gamma;
  
  if(bInteraction)
    gamma = as<NumericMatrix>(results["gamma"]);
  
  for(int row = 0; row < nout - 1; row++){
    for(int loc = 0; loc < nloc; loc++){
      for(int pop = 0; pop < npop; pop++){
        double w = alpha(row, loc) + beta(row, pop) + (bInteraction ? gamma(row, loc * npop + pop) : 0);
        fst(row, loc * npop + pop) = std::exp(w) / (1 + std::exp(w));
      }
    }
  }
  
  return fst;
}
    
