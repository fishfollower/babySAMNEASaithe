#include <TMB.hpp>
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}

template <class Type>
Type itrans(Type x){ // scaled ilogit
  return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);
}

template<class Type>
vector<Type> ssbFUN(matrix<Type> logN, matrix<Type> logFF, matrix<Type> M, matrix<Type> SW, matrix<Type> MO, matrix<Type> PF, matrix<Type> PM){
  int nrow = logN.rows();
  int ncol = logN.cols();
  vector<Type> ret(nrow);
  ret.setZero();
  for(int y=0; y<nrow; ++y){
    for(int a=0; a<ncol; ++a){
      ret(y) += SW(y,a)*MO(y,a)*exp(logN(y,a))*exp(-PF(y,a)*exp(logFF(y,a))-PM(y,a)*M(y,a));
    }
  }
  return ret;
}

template<class Type>
vector<Type> xExpA(vector<Type> x, Eigen::SparseMatrix<Type> A, int Nstep=20, bool uniformization=false){
  Type rho=0;
  if(uniformization){
    vector<Type> d = A.diagonal();
    using TMBad::min;
    for(int i=0;i<d.size();++i)rho=min(rho,d(i));
    Eigen::SparseMatrix<Type> D(A.rows(),A.cols());
    for(int i=0; i<D.rows();++i){D.coeffRef(i,i)=-rho;}
    A = A+D;
  }
  matrix<Type> Yrow(1,x.size()), lastTerm(1,x.size()); 
  Yrow.row(0)=x;
  lastTerm.row(0)=x;
  for(int i=1; i<=Nstep; ++i){
    lastTerm=lastTerm*A/i;
    Yrow+=lastTerm;
  }
  vector<Type> y=Yrow.row(0);
  if(uniformization){
    y=exp(rho)*y;
  }
  return y;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs)  
  DATA_IMATRIX(aux)
  DATA_IMATRIX(idx1)
  DATA_IMATRIX(idx2)
  DATA_INTEGER(minYear)
  DATA_INTEGER(minAge)    
  DATA_IVECTOR(fleetTypes)
  DATA_VECTOR(sampleTimes)     
  DATA_VECTOR(year)
  DATA_VECTOR(age)
  DATA_MATRIX(M)
  DATA_MATRIX(SW)
  DATA_MATRIX(MO)
  DATA_MATRIX(PF)
  DATA_MATRIX(PM)

  DATA_INTEGER(srmode)
  DATA_INTEGER(fcormode)
  DATA_IVECTOR(keyF)       
  DATA_IMATRIX(keyQ)   
  DATA_IMATRIX(keySd)
  DATA_IVECTOR(covType)
  DATA_IMATRIX(keyIGAR)
  DATA_IVECTOR(noParUS)
  DATA_INTEGER(mode)
  int nobs=obs.size();    
  int nrow=M.rows();
  int ncol=M.cols();
  
  PARAMETER(logsdR)
  PARAMETER(logsdS)
  PARAMETER_VECTOR(logsdF)
  PARAMETER_VECTOR(rickerpar)
  PARAMETER_VECTOR(transRhoF)  
  PARAMETER_VECTOR(bhpar)    
  PARAMETER_VECTOR(logQ)
  PARAMETER_VECTOR(logsd)
  PARAMETER_VECTOR(logIGARdist)
  PARAMETER_VECTOR(parUS)
  PARAMETER_MATRIX(logN)
  PARAMETER_MATRIX(logF)    
  PARAMETER_VECTOR(missing)
  Type sdR = exp(logsdR);
  Type sdS = exp(logsdS);
  vector<Type> sdF = exp(logsdF);      
  vector<Type> sd=exp(logsd);
  Type rhoF=0;
  
  //patch missing 
  int idxmis=0; 
  for(int i=0; i<nobs; i++){
    if(isNA(obs(i))){
      obs(i)=exp(missing(idxmis++));
    }    
  }

  matrix<Type> logFF(nrow,ncol);
  for(int i=0; i<nrow; ++i){
    for(int j=0; j<ncol; ++j){
      logFF(i,j)=logF(i,keyF(j));
    }
  }

  int Adim = 3*ncol;  
  matrix<Type> firstState(nrow,Adim);
  firstState.setZero();
  matrix<Type> lastState(nrow,Adim);  
  lastState.setZero();
  
  // non-sparse ed ////////////////////////////////////////////////
  if(mode==0){
    vector< matrix<Type> > Avec(nrow);
     
    for(int i=0; i<nrow; ++i){
      matrix<Type> AA(Adim,Adim); AA.setZero();
      for(int j=0; j<ncol; ++j){
        AA(j,j+ncol) = exp(logFF(i,j));
        AA(j,j+2*ncol) = M(i,j);
        AA(j,j) = -exp(logFF(i,j))-M(i,j); // negative row sum
      }
      Avec(i) = AA;
    }
    

    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
        firstState(i,j)=exp(logN(i,j));
      }
      lastState.row(i)=firstState.row(i)*atomic::expm(Avec(i)); 
    }
    REPORT(firstState);
    REPORT(lastState);
  }
  //  std::cout<<lastState.row(5)<<std::endl;

  
  //  sparse ed  ///////////////////////////////////////////////////
  if(mode==1){
    vector< Eigen::SparseMatrix<Type> > Avec(nrow);
    
    for(int i=0; i<nrow; ++i){
      Eigen::SparseMatrix<Type> AA(Adim,Adim);
      for(int j=0; j<ncol; ++j){
        AA.coeffRef(j,j+ncol) = exp(logFF(i,j));
        AA.coeffRef(j,j+2*ncol) = M(i,j);
        AA.coeffRef(j,j) = -exp(logFF(i,j))-M(i,j); // negative row sum
      }
      // to cheat TMB (expm seems not to like sparse-zero diag elements)
        for(int j=ncol; j<AA.cols(); ++j){
          AA.coeffRef(j,j) = 0; 
        }
        Avec(i) = AA;
      }
  

    sparse_matrix_exponential::config<Type> myconfig = sparse_matrix_exponential::config<Type>();
    myconfig.Nmax = 1000;
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
        firstState(i,j)=exp(logN(i,j));
      }
      Eigen::SparseMatrix<Type> tmp=Avec(i);
      sparse_matrix_exponential::expm_generator<Type> eg(tmp,myconfig);
      vector<Type> t2=firstState.row(i);
      lastState.row(i)=eg(t2).array();
    }
  }
  //  //std::cout<<lastState.row(5)<<std::endl;

  // sparse 2. ed  /////////////////////////////////////////////////// 
  //  vector< Eigen::SparseMatrix<Type> > Avec(nrow);
  //  
  //  int Adim = 3*ncol;
  //  for(int i=0; i<nrow; ++i){
  //    Eigen::SparseMatrix<Type> AA(Adim,Adim);
  //    for(int j=0; j<ncol; ++j){
  //      AA.coeffRef(j,j+ncol) = exp(logFF(i,j));
  //      AA.coeffRef(j,j+2*ncol) = M(i,j);
  //      AA.coeffRef(j,j) = -exp(logFF(i,j))-M(i,j); // negative row sum
  //    }
  //    Avec(i) = AA;
  //  }
  //  
  //  matrix<Type> firstState(nrow,Adim);
  //  firstState.setZero();
  //  matrix<Type> lastState(nrow,Adim);  
  //  lastState.setZero();
  //  int Nstep = 20;
  //  for(int i=0; i<nrow; ++i){
  //    for(int j=0; j<ncol; ++j){
  //      firstState(i,j)=exp(logN(i,j));
  //    }
  //    lastState.row(i)=firstState.row(i);
  //    matrix<Type> lastTerm=lastState.row(i);
  //    Eigen::SparseMatrix<Type> tmp=Avec(i);
  //    for(int ii=1; ii<=Nstep; ++ii){
  //      lastTerm=lastTerm*tmp/ii;
  //      lastState.row(i)+=lastTerm;
  //    }
  //  }


  // sparse 3. ed  ///////////////////////////////////////////////////
  if(mode==2){
    vector< Eigen::SparseMatrix<Type> > Avec(nrow);
    
    for(int i=0; i<nrow; ++i){
      Eigen::SparseMatrix<Type> AA(Adim,Adim);
      for(int j=0; j<ncol; ++j){
        AA.coeffRef(j,j+ncol) = exp(logFF(i,j));
        AA.coeffRef(j,j+2*ncol) = M(i,j);
        AA.coeffRef(j,j) = -exp(logFF(i,j))-M(i,j); // negative row sum
      }
      Avec(i) = AA;
    }
    
    for(int i=0; i<nrow; ++i){
      for(int j=0; j<ncol; ++j){
        firstState(i,j)=exp(logN(i,j));
      }
      lastState.row(i)=xExpA(vector<Type>(firstState.row(i)),Avec(i)); 
    }
  }
  Type jnll=0;
  
  vector<Type> ssb = ssbFUN(logN,logFF,M,SW,MO,PF,PM);

  // ############################################# N part
  Type predN, thisSSB;
  for(int y=1; y<nrow; ++y){

    if((y-minAge)>=0){
      thisSSB=ssb(y-minAge);
    }else{
      thisSSB=ssb(0); // use first in beginning       
    } 

    switch(srmode){
      case 0:
        predN = logN(y-1,0);
      break;

      case 1:
	predN = rickerpar(0)+log(thisSSB)-exp(rickerpar(1))*thisSSB;
      break;

      case 2:
        predN = bhpar(0)+log(thisSSB)-log(1.0+exp(bhpar(1))*thisSSB);
      break;

      default:
	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
	exit(EXIT_FAILURE);
      break;
    }
    jnll += -dnorm(logN(y,0),predN,sdR,true);
  }
  
  for(int y=1; y<nrow; ++y){
    for(int a=1; a<ncol; ++a){
      predN=log(lastState(y-1,a-1));
      if(a==(ncol-1)){
        predN=log(lastState(y-1,a-1)+lastState(y-1,a));
      }
      jnll += -dnorm(logN(y,a),predN,sdS,true);
    }
  }

  // ############################################# F part
  matrix<Type> SigmaF(logF.cols(),logF.cols());
  SigmaF.setZero();

  switch(fcormode){

    case 0:
      SigmaF.diagonal() = sdF*sdF;
    break;

    case 1:
      SigmaF.diagonal() = sdF*sdF;
      rhoF = itrans(transRhoF(0));
      for(int i=0; i<logF.cols(); ++i){
        for(int j=0; j<i; ++j){
          SigmaF(i,j) = rhoF*sdF(i)*sdF(j);
          SigmaF(j,i) = SigmaF(i,j); 
        }
      }
    break;

    case 2:
      SigmaF.diagonal() = sdF*sdF;
      rhoF = itrans(transRhoF(0));
      for(int i=0; i<logF.cols(); ++i){
        for(int j=0; j<i; ++j){
          SigmaF(i,j) = pow(rhoF,Type(i-j))*sdF(i)*sdF(j);
          SigmaF(j,i) = SigmaF(i,j); 
        }
      }
    break;

    default:
      std::cout<<"Error: This cormode not implemented yet."<<std::endl;
      exit(EXIT_FAILURE);
    break;
  }
  
  density::MVNORM_t<Type> nldens(SigmaF);
  for(int y=1; y<nrow; ++y){
    jnll += nldens(logF.row(y)-logF.row(y-1));
  }
  /// obs part
  
  vector<Type> logPred(nobs);
  Type Z;
  int y, a, f;
  for(int i=0; i<nobs; ++i){
    y=aux(i,0)-minYear;
    f=aux(i,1)-1;
    a=aux(i,2)-minAge;
    Z=exp(logFF(y,a))+M(y,a);
    switch(fleetTypes(f)){
      case 0:
        logPred(i) = log(lastState(y,a+ncol));//logN(y,a)-log(Z)+log(1-exp(-Z))+logFF(y,a);
      break;
      
      case 2:
        logPred(i) = logQ(keyQ(f,a))+logN(y,a)-Z*sampleTimes(f); //logQ(keyQ(f,a))+log(xExpA(vector<Type>(firstState.row(y)),  Eigen::SparseMatrix<Type>(Avec(y)*sampleTimes(f))))(a); // 
      break;

      default:
        std::cout<<"Error: This fleet type not implemented yet."<<std::endl;
        exit(EXIT_FAILURE);
      break;
    }
  }

  vector< density::MVNORM_t<Type> >  nllVec(fleetTypes.size());
  vector< matrix<Type> >  Svec(fleetTypes.size());

  for(int f=0; f<idx1.rows(); ++f){
    int d=0;
    for(int a=0; a<keySd.cols(); ++a){
      if(!isNAINT(keySd(f,a)))d++;
    };
    int thisdim=d;
    matrix<Type> S(thisdim,thisdim);
    S.setZero();
    if(covType(f)==0){
      d=0;
      for(int a=0; a<keySd.cols(); ++a){
        if(!isNAINT(keySd(f,a))){
          S(d,d)=sd(keySd(f,a))*sd(keySd(f,a));
          d++;
        }
      };
    }
    if(covType(f)==1){
      vector<Type> dist(thisdim);
      dist.setZero();
      d=1;
      for(int a=0; a<keyIGAR.cols(); ++a){
         if(!isNAINT(keyIGAR(f,a))){
           dist(d)=dist(d-1)+exp(logIGARdist(keyIGAR(f,a)));
           d++;
         }
      }
      vector<Type> sdvec(thisdim);
      d=0;
      for(int a=0; a<keySd.cols(); ++a){
        if(!isNAINT(keySd(f,a))){
          sdvec(d++)=sd(keySd(f,a));
        }
      }
      for(int i=0; i<S.rows(); ++i){
        for(int j=0; j<i; ++j){
          S(i,j) = pow(0.5,dist(i)-dist(j))*sdvec(i)*sdvec(j);
          S(j,i) = S(i,j);
        }
        S(i,i)=sdvec(i)*sdvec(i);
      }
    }
    if(covType(f)==2){
      vector<Type> sdvec(thisdim);
      d=0;
      for(int a=0; a<keySd.cols(); ++a){
        if(!isNAINT(keySd(f,a))){
          sdvec(d++)=sd(keySd(f,a));
        }
      }
      int from=0;
      for(int i=0; i<f; ++i) from+=noParUS(i);
      vector<Type> thispar=parUS.segment(from,noParUS(f));
      density::UNSTRUCTURED_CORR_t<Type> USDENS(thispar);
      matrix<Type> DD(S.rows(),S.cols());
      DD.setZero();
      for(int i=0; i<DD.rows(); ++i)DD(i,i)=sdvec(i);
      S = DD*(USDENS.cov())*DD;
    }

    nllVec(f).setSigma(S);
    Svec(f)=S;
  }
  
  for(int f=0; f<idx1.rows(); ++f){
    for(int y=0; y<idx1.cols(); ++y){
      if(!isNAINT(idx1(f,y))){
        vector<Type> o=obs.segment(idx1(f,y),idx2(f,y)-idx1(f,y)+1);
        vector<Type> p=logPred.segment(idx1(f,y),idx2(f,y)-idx1(f,y)+1);
        jnll += nllVec(f)(log(o)-p);
      }
    }
  }
  vector<Type> logObs=log(obs);
  REPORT(logObs)
  REPORT(logPred)
  REPORT(Svec);
  ADREPORT(ssb);
  return jnll;
}
