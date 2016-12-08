#include <algorithm>
#include <random>
#include <Rcpp.h>

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
    
    // input variables which were constants and so may be rarely changed
    int m_nNumOut; // default 2001
    
    // Random number generator stuff
    mt19937 m_mtEngine;
    normal_distribution<double> m_rngZ;
    uniform_real_distribution<double> m_rngU;
    
    // internal variables
    
    vector<double> alpha, beta, gamma;
    vector<double> lp0, lp1;
    vector< vector<double> > pi;
    bool bIllegal;
    vector<double> m_vdFst;
    vector<int> m_vnNumZeros;
    double m_dMeanFst;

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
    
    void initHyperparams(){
        alpha.resize(m_nLoci);
        
        alpha[0] = rnorm(alphaMu, alphaSigma);
        for(int i = 1; i < m_nLoci; i++)
            alpha[i] = (cor * alpha[i-1]) + rnorm() * srk * alphaSigma;
        
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
                           vector<double>& lp, vector<double>& pa){
        double lpd = 0;
        
        for(int j = 0; j < m_nPops; j++)
            lpd += logdnorm(b[j] , betaMu , betaSigma);     /* PRIOR FOR b_j */
        
        for(int i = 0;i < m_nLoci; i++){
            if(i==0){
                lp[i] = logdnorm(a[i], alphaMu, alphaSigma);	/* PRIOR FOR a_i */
            }else if(i > 0 && i < m_nLoci){
                lp[i]= logdnorm(a[i], cor * a[i-1], srk * alphaSigma);
                pa[i]= lp[i];
            }

            pa[m_nLoci] = 0;

            for (int j = 0; j < m_nPops; j++){
                if (m_popSums[i][j] > 0){
                    lp[i] += logMultinomDirichlet( calcFst(a[i] + b[j]), i, j);   /* LIKELIHOOD */
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
    
    double logDirichlet(int k, vector<double>& par, vector<double> vec){
        double sum;
        double konst;
        
        sum = 0.0;
        konst = 0.0;
        
        for(int j = 0; j < k; j++){
            sum += (par[j] - 1.0) * std::log(vec[j]);
            sum -= std::lgamma(par[j]);
            konst += par[j];
        }
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
    
    void pargen(int i, const vector<double>& oldpi){
        vector<double> fpar, rpar, tvec;
        int k;
        double cum = 0.0;
        int nA = m_numAlleles[i];
        
        fpar.resize(nA);
        rpar.resize(nA);
        tvec.resize(nA);
        
        for(k = 0; k < nA; k++){
            tvec[k] = oldpi[k];
            
            if(tvec[k] < 1.0e-3){
                tvec[k] = 1.0e-3;
            }
            
            fpar[k] =  pSigma * tvec[k];
            pi[i][k] = rgamma(fpar[k], 1.0)[0]; // rgamma is from Rcpp so it is vectorized
            
            if(pi[i][k] == 0.0){
                bIllegal = 1;
                return;
            }
            
            cum += pi[i][k];
        }
        
        for(k = 0; k < nA; k++){
            pi[i][k] /= cum;
            tvec[k] = pi[i][k];
            
            if(tvec[k] < 1.0e-3){
                tvec[k] = 1.0e-3;
            }
            rpar[k] = pSigma * tvec[k];
        }
        
        /* forward and reverse  log-Hastings terms */
        ftp = logDirichlet(nA, fpar, pi[i]);
        rtp = logDirichlet(nA, rpar, oldpi);
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
      
      for(int i = 0; i < m_nLoci; i++){
        for(int k = 0; k < m_numAlleles[i]; k++){
          Rprintf("%6.3f", pi[i][k]);
        }
        Rprintf("\n");
      }
      Rprintf("\n");
    }
    
    
    double quantileFst(double p){
        double idx = p * m_nNumOut + 0.5;
        if(idx - std::floor(idx) < 0.1){
            return m_vdFst[(int)idx];
        } //else{
        int k = (int)idx;
        return m_vdFst[k] + fmod(idx, 1.0) * (m_vdFst[k + 1] - m_vdFst[k]);
    }
    
    double runif(void){
        return m_rngU(m_mtEngine);
    }
    
    double rnorm(double mu = 0, double sigma = 1){
        return sigma * m_rngZ(m_mtEngine) + mu;
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
        
        double usa = alphaSigma * uSigma;
        double usg = gammaSigma * uSigma;
        
        for (i = 0; i < m_nLoci; i++){
            vector<double>  oldpi = pi[i];
            
            pargen(i, oldpi);
            
            alphaDash[i] = rnorm(a[i], usa);
            
            if(i == 0){
                lp1[0] = logdnorm(alphaDash[0] , alphaMu , alphaSigma);
                pa1[0] = logdnorm(alphaDash[0] , alphaMu , alphaSigma);
                pa1[1] = logdnorm(a[1] , cor * alphaDash[0] , srk * alphaSigma);
            }
            else if(i > 0 &&i < m_nLoci - 1){
                lp1[i]= logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i]= logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i + 1]= logdnorm(a[i + 1] , cor * alphaDash[i] , srk * alphaSigma);
            }
            else if(i == m_nLoci - 1){
                lp1[i]= logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
                pa1[i]= logdnorm(alphaDash[i] , cor * a[i - 1] , srk * alphaSigma);
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
        
        double usb = betaSigma * uSigma / std::sqrt(m_nLoci);
        
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
        
        /* CORRELATION PARAMETER (FIXED) BETWEEN ADJACENT MARKERS */
        
        cor = 0; // the value of cor determines the shrinkage value
        // associated with the sd of each a_i so that the
        // overall - expected - sd is alphaSigma (defined above).
        // The shrinkage value (srk) is defined below
        
        srk = 1;
        
        // Interaction is switched OFF by default
        m_bGammaSwitch = false;
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
    
    void printRunInfo(void){
        Rprintf("****** Settings for MCMC algorithm ******\n\n");
        Rprintf("Burn-in length (iterations): %ld\n", m_nDiscard);
        Rprintf("Thinning interval: %ld\n", m_nKeep);
        Rprintf("Number of iterations (post burn-in): %ld\n" , m_nNumIt);
        Rprintf("Number of outputs: %ld\n", m_nNumOut);
        
        Rprintf("Using Normal prior for hyperparameters a_i with mean %.2f, sd %.2f\n", alphaMu, alphaSigma);
        Rprintf("Using Normal prior for hyperparameters b_j with mean %.2f, sd %.2f\n", betaMu, betaSigma);
        
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
    }
    
    void setData(List data){
        m_nPops = data["nPops"];
        m_nLoci = data["nLoci"];
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
            
            Rprintf("Using Normal prior for hyperparameters a_i with mean %.2f, sd %.2f\n",alphaMu,alphaSigma);
            Rprintf("Using Normal prior for hyperparameters b_j with mean %.2f, sd %.2f\n",betaMu,betaSigma);
            
            if (m_bGammaSwitch){
                Rprintf("Using interaction terms in model\n");
                Rprintf("Using Normal prior for hyperparameters g_ij with mean %.2f, sd %.2f\n\n",gammaMu,gammaSigma);
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
        
        // Initialise allele frequencies
        initAlleleFreqs();
        //printHyperparams();
        
        /* INITIALISE LOG-POSTERIOR DENSITY VALUE (z) */
        
        vector<double> lp, pa, pg;
        lp.resize(m_nLoci);
        pa.resize(m_nLoci + 1);

        if(m_bGammaSwitch)
            pg.resize(m_nLoci);

        double logPosteriorDensity = initLogPostDens(alpha, beta, gamma, lp, pa);


        /* RUN METROPOLIS ALGORITHM (for numit+discard iterations) */

        int jmp1 = 0;
        int jmp2 = 0;
        
        NumericMatrix postAlpha(m_nNumOut, m_nLoci);
        NumericMatrix postBeta(m_nNumOut, m_nPops);
        NumericMatrix postGamma(m_nNumOut, m_nLoci * m_nPops);
        NumericVector lpd(m_nNumOut);
        
        int sumAlleles = std::accumulate(m_numAlleles.begin(), m_numAlleles.end(), 0.0);
        NumericMatrix postP(m_nNumOut, sumAlleles);
        
        int ctr = 0;
        Progress p1(m_nDiscard, true);
         
        for (int nCurrentIteration = -m_nDiscard; nCurrentIteration < m_nNumIt; nCurrentIteration++){
            bIllegal = false;
            
            jmp1 += update_beta(alpha, beta, gamma, logPosteriorDensity, lp, pa);
            jmp2 += update_alpha(alpha, beta, gamma, logPosteriorDensity, lp, pa, pg);
            
            if(nCurrentIteration < 0)
              p1.increment();
        
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
                if (((nCurrentIteration + m_nDiscard) > 0) && ((nCurrentIteration+m_nDiscard)/m_nAcceptanceRateGap)*m_nAcceptanceRateGap==(nCurrentIteration+m_nDiscard)){
                  double percentDone = 100.0 * nCurrentIteration / (double) m_nNumIt; 
                    Rprintf("Percent: %5.2f Iter: %8d, Accpt. rate: %12.6f%12.6f\n", percentDone, nCurrentIteration, ((double)(jmp1)/m_nAcceptanceRateGap), ((double)(jmp2)/m_nAcceptanceRateGap/m_nLoci));
                    jmp1=0; jmp2=0;
                }
            }
        }
        
        List results;
        
        results["lpd"] = lpd;
        results["alpha"] = postAlpha;
        results["beta"] = postBeta;
        results["p"] = postP;
        if(m_bGammaSwitch)
          results["gamma"] = postGamma;
        
        return results;
    }
};
    
    
    
    RCPP_MODULE(BayesFst) {
        using namespace Rcpp;
        
        class_<BayesFst>( "BayesFst")
        .default_constructor("Standard constructor")
        .method("printData", &BayesFst::printData)
        .method("printFstSummary", &BayesFst::printFstSummary)
        .method("run", &BayesFst::run)
        .method("setData", &BayesFst::setData)
        .method("setPriorParameters", &BayesFst::setPriorParameters)
        .method("setRunParameters", &BayesFst::setRunParameters)
        .property("interaction", &BayesFst::getInteraction, &BayesFst::setInteraction)
        ;
    }
    
    