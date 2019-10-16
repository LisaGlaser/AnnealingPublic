#include "TwoSphere.h"

void TwoSphere::setsize(int mysize){
  sizeS2 = mysize;
  spinorsize = 2*sizeS2;
  totalsize=2*spinorsize;
  eigenSpinors = produceEigenSpinors();
  PBott=producePMat();
  gamma = produceGamma();
  DiracS2 = produceDiracS2();
  PosDiracS2 = producePosDiracS2();
}

int TwoSphere::getsize(){
  return sizeS2;
}

// convention: l,lp are odd integers (twice the l-values), m,mp are twice the corresponding m-values, false means positive and true means negative eigenvals of D (so that the {-1,1}-valued sign is (-1)^sgn)

double TwoSphere::getSU2GeneratorACoefficients(int dl, int dm, bool sgn, int dlp, int dmp, bool sgnp){
// we return the coefficients <lp,mp|_{sgnp} a |l,m>_{sgn}
// numbers from Dabrowski et al 2005, but our |l,m>_+ is their |l,m>_+ + |l,m>-- and
// our |l,m>_- is their - |l,m>_+ + |l,m>_-.
  auto dabrowskicoeff = [](int dl, int dm, bool sgn, int dlp, int dmp, bool sgnp)->double{
    if ( dlp==dl && dmp == dm+2 && sgn==sgnp ){
      return 2*pow(-1,sgn)*sqrt((dl+dm+2)*(dl-dm)/4)/(dl*(dl+2));
    }
    else if ( dlp == dl+2 && dmp == dm+2 && sgn==sgnp ){
      return sqrt((dl+dm+2)*(dl+dm+4)/4)/(dl+2);
    }
    else if ( dlp == dl -2 && dmp == dm+2 && sgn== sgnp ){
      return -1*sqrt((dl-dm)*(dl-dm-2)/4)/(dl);
    }
    else{
      return 0;
    }
  };
  return (pow(-1,int(sgn+sgnp))*dabrowskicoeff(dl,dm,false,dlp,dmp,false) + dabrowskicoeff(dl,dm,true,dlp,dmp,true))/2;
}

double TwoSphere::getSU2GeneratorBCoefficients(int dl, int dm, bool sgn, int dlp, int dmp, bool sgnp){
// we return the coefficients <lp,mp|_{sgnp} b |l,m>_{sgn}
// numbers from Dabrowski et al 2005
  auto dabrowskicoeff = [](int dl, int dm, bool sgn, int dlp, int dmp, bool sgnp)->double{
    if ( dlp==dl && dmp == dm && sgn==sgnp ){
      return pow(-1,sgn)*((dl-dm+2)*(dl+dm)/4 - (dl-dm)*(dl+dm+2)/4)/(dl*(dl+2));
    }
    else if ( dlp == dl+2 && dmp == dm && sgn==sgnp ){
      return -1*sqrt((dl-dm+2)*(dl+dm+2)/4)/(dl+2);
    }
    else if ( dlp == dl -2 && dmp == dm && sgn== sgnp ){
      return -1*sqrt((dl-dm)*(dl+dm)/4)/(dl);
    }
    else{
      return 0;
    }
  };
  return (pow(-1,int(sgn+sgnp))*dabrowskicoeff(dl,dm,false,dlp,dmp,false) + dabrowskicoeff(dl,dm,true,dlp,dmp,true))/2;
}

s2mode_vect TwoSphere::produceEigenSpinors(){
  s2mode_vect eS;
  eS.resize(0);
  s2mode eigenvector;
  int dl,dm;
  for ( dl = 1; eS.size()<spinorsize; dl=dl+2 ){
    // std::cout<<"dl is "<< dl << std::endl;
    for ( dm = dl; dm >= -dl; dm = dm-2 ){
      // std::cout<<"dm is "<< dm << std::endl;
      for ( bool sgn : { false, true }){
        // std::cout<<"sgn is "<< sgn << std::endl;
        eigenvector = std::make_tuple(dl,dm,sgn);
        eS.push_back(eigenvector);
        }
    }
  }
  sort(eS.begin(), eS.end(), TwoSphere::eigenSort);
  eS.resize(spinorsize);
  return eS;
}

bool TwoSphere::eigenSort(s2mode i,s2mode j){
  // minus before plus,
  if ( std::get<2>(i) != std::get<2>(j) ){
    return std::get<2>(i);
  }
  // then ascending in l,
  else if ( std::get<0>(i) != std::get<0>(j)){
    return std::get<0>(i) < std::get<0>(j);
  }
  // then ascending in m
  else if ( std::get<1>(i) != std::get<1>(j)){
    return std::get<1>(i) < std::get<1>(j);
  }
  else {
    return false;
  }
}


Eigen::MatrixXcd TwoSphere::setSpinorsizedMatrix(double (coeff)(int,int,bool,int,int,bool)){
  int l,m,lp,mp;
  bool sgn,sgnp;
  Eigen::MatrixXcd mat(spinorsize,spinorsize);
  for ( int i = 0; i < eigenSpinors.size(); i++){
    for ( int j = 0; j < eigenSpinors.size(); j++){
      l = std::get<0>(eigenSpinors[i]);
      m = std::get<1>(eigenSpinors[i]);
      sgn=std::get<2>(eigenSpinors[i]);
      lp = std::get<0>(eigenSpinors[j]);
      mp = std::get<1>(eigenSpinors[j]);
      sgnp=std::get<2>(eigenSpinors[j]);
      mat(i,j) = coeff(l,m,sgn,lp,mp,sgnp);
    }
  }
  return mat;
}


Eigen::MatrixXcd TwoSphere::produceLeftUpperPMat(){
  auto coeff = [](int l,int m,bool sgn,int lp,int mp ,bool sgnp)->double {
    return (int(l==lp && m == mp && sgn == sgnp) - getSU2GeneratorBCoefficients(l,m,sgn,lp,mp,sgnp))/2;
  };
  Eigen::MatrixXcd mat = setSpinorsizedMatrix(coeff);
  return mat;
}

Eigen::MatrixXcd TwoSphere::produceRightUpperPMat(){
  auto coeff = [](int l,int m,bool sgn,int lp,int mp ,bool sgnp)->double {
    return getSU2GeneratorACoefficients(lp,mp,sgnp,l,m,sgn)/2;
  };
  Eigen::MatrixXcd mat = setSpinorsizedMatrix(coeff);
  return mat;
}

Eigen::MatrixXcd TwoSphere::produceLeftLowerPMat(){
  auto coeff = [](int l,int m,bool sgn,int lp,int mp ,bool sgnp)->double {
    return getSU2GeneratorACoefficients(l,m,sgn,lp,mp,sgnp)/2;
  };
  Eigen::MatrixXcd mat = setSpinorsizedMatrix(coeff);
  return mat;
}

Eigen::MatrixXcd TwoSphere::produceRightLowerPMat(){
  auto coeff = [](int l,int m,bool sgn,int lp,int mp ,bool sgnp)->double {
    return (int(l==lp && m == mp && sgn == sgnp) + getSU2GeneratorBCoefficients(l,m,sgn,lp,mp,sgnp))/2;
  };
  Eigen::MatrixXcd mat = setSpinorsizedMatrix(coeff);
  return mat;
}

Eigen::MatrixXcd TwoSphere::producePMat(){
  Eigen::MatrixXcd mat(totalsize,totalsize);
  mat.topLeftCorner(spinorsize,spinorsize) = produceLeftUpperPMat();
  mat.bottomRightCorner(spinorsize,spinorsize) = produceRightLowerPMat();
  mat.topRightCorner(spinorsize,spinorsize) = produceRightUpperPMat();
  mat.bottomLeftCorner(spinorsize,spinorsize) = produceLeftLowerPMat();
  return mat;
}

Eigen::MatrixXcd TwoSphere::RelativeTrace(Eigen::MatrixXcd mat){
  Eigen::MatrixXcd trmat(spinorsize,spinorsize);
  if ( mat.rows() != totalsize || mat.cols() != totalsize){
    std::cout << "Wrong matrix size in relativetrace!" << std::endl;
  }
  trmat = 2*(mat.topLeftCorner(spinorsize,spinorsize)+mat.bottomRightCorner(spinorsize,spinorsize));
  return trmat;
}

Eigen::MatrixXcd TwoSphere::doubleMatrix(Eigen::MatrixXcd mat){
  Eigen::MatrixXcd doublemat;
  doublemat = Eigen::MatrixXcd::Zero(totalsize,totalsize);
  if ( mat.rows() != spinorsize || mat.cols() != spinorsize ){
    std::cout << "Wrong matrix sizes in doubleMatrix!" << std::endl;
  }
  doublemat.topLeftCorner(spinorsize,spinorsize) = mat;
  doublemat.bottomRightCorner(spinorsize,spinorsize) = mat;
  return doublemat;
}

Eigen::MatrixXcd TwoSphere::producePosDiracS2(){
  Eigen::MatrixXcd dirac = produceDiracS2();
  return dirac.bottomRightCorner(sizeS2,sizeS2);
}
Eigen::MatrixXcd TwoSphere::getPosDiracS2sqrt(){

  Eigen::MatrixXcd dirac = produceDiracS2();
  Eigen::MatrixXcd sqrtD=dirac.bottomRightCorner(sizeS2,sizeS2);
  for (int j=0;j<sqrtD.cols();j++)
    {
        for (int i=0;i<sqrtD.rows();i++)
        {

            sqrtD(j,i)=sqrt(sqrtD(j,i));

        }
    }
  return sqrtD;
}


Eigen::MatrixXcd TwoSphere::produceDiracS2(){
  auto coeff = [](int l, int m, bool sgn, int lp, int mp, bool sgnp)->double{
    // D|l,m>_\pm = (l+1/2)|l,m>_{\mp}
    auto dabrowskicoeff = [](int l, int m, bool sgn, int lp, int mp, bool sgnp)->double{
      return int(l==lp && m == mp && sgn != sgnp)*(l+1)/2;
    };
    // translate bases as before
    return (pow(-1,int(sgn+sgnp))*dabrowskicoeff(l,m,false,lp,mp,false) + pow(-1,int(sgnp))*dabrowskicoeff(l,m,true,lp,mp,false) + pow(-1,int(sgn))*dabrowskicoeff(l,m,false,lp,mp,true) + dabrowskicoeff(l,m,true,lp,mp,true))/2;
  };
  Eigen::MatrixXcd dirac = setSpinorsizedMatrix(coeff);
  return dirac;
}

Eigen::MatrixXcd TwoSphere::oneSidedEquationLHS(Eigen::MatrixXcd D, Eigen::MatrixXcd p){
  if ( D.rows() != spinorsize || D.cols() != spinorsize || p.rows() != totalsize || p.cols() != totalsize ){
    std::cout << "Wrong matrix sizes in oneSidedEquationLHS!" << std::endl;
    return Eigen::MatrixXcd::Zero(totalsize,totalsize);
  }
  Eigen::MatrixXcd Dp = doubleMatrix(D);
  Eigen::MatrixXcd one = Eigen::MatrixXcd::Identity(totalsize,totalsize);
  Eigen::MatrixXcd commutator = Dp*p - p*Dp;
  Eigen::MatrixXcd beforetrace=(p-one/2)*commutator*commutator;
  return RelativeTrace(beforetrace);
}

Eigen::MatrixXcd TwoSphere::produceGamma(){
  Eigen::MatrixXcd gammatmp;
  gammatmp = Eigen::MatrixXcd::Zero(spinorsize,spinorsize);
  gammatmp.topRightCorner(sizeS2,sizeS2) = - Eigen::MatrixXcd::Identity(sizeS2,sizeS2);
  gammatmp.bottomLeftCorner(sizeS2,sizeS2) = - Eigen::MatrixXcd::Identity(sizeS2,sizeS2);
  return gammatmp;
}

Eigen::MatrixXcd TwoSphere::getGamma(){
  return gamma;
}

Eigen::MatrixXcd TwoSphere::getPMat(){
  return PBott;
}

s2mode_vect TwoSphere::getEigenSpinors(){
  return eigenSpinors;
}

Eigen::MatrixXcd TwoSphere::getDiracS2(){
  return DiracS2;
}

Eigen::MatrixXcd TwoSphere::getPosDiracS2(){
  return PosDiracS2;
}

double TwoSphere::oneSidedEquationFrobSq(Eigen::MatrixXcd D, Eigen::MatrixXcd p){
  if ( D.rows() != spinorsize || D.cols() != spinorsize || p.rows() != totalsize || p.cols() != totalsize ){
    std::cout << "Wrong matrix sizes in oneSidedEquationLHS!" << std::endl;
    return 0;
  }
  Eigen::MatrixXcd diff = oneSidedEquationLHS(D,p)-gamma;
  double frob = diff.squaredNorm();
  return frob;
}
