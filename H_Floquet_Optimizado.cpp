//--------------libraries----------------------------
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <matrixtypes.h> 
#include <lapack.h>
#include <iostream> 
#include <fstream> 
#include <iomanip>

//--------------------------------------------------
using namespace boost::numeric::ublas;
using namespace std;
using namespace ula;
using boost::math::quadrature::trapezoidal;
using namespace  boost::math::constants;

//--------------------One Particle Bose-Hubbard Hamiltonian----------------
class Hamiltonian{
  public:
    int n_sites;
    double mu, J, q, T;
    Complex iu;
    ComplexMatrix Identity;

    Hamiltonian(int _n_sites, double _mu, double _J, double _q, double _T){
      n_sites=_n_sites;
      mu=_mu; J=_J; q=_q; T=_T;
      iu=complex<double>(0,1);
      identity_matrix<double> id(4*n_sites);
      Identity=id;
    }
      
    ComplexMatrix FreeHamiltonian(){
      int sizeH=4*n_sites;  // # de estados y por tanto dimensión de la matriz

      RealMatrix FH(sizeH,sizeH,0.0);  

      // Términos diagonales (Chemical potential)

      for(int j=0; j<sizeH ; j=j+4){

        FH(j+1,j+1)=-mu;
        FH(j+2,j+2)=-mu;
        FH(j+3,j+3)=-mu;

      }
        
      // Hopping 

      if(n_sites>2){
        for(int i=0; i<sizeH-4 ; i=i+4){
          FH(i+1,i+5) = -J;  FH(i+5,i+1) = -J;
          FH(i+2,i+6) = -J;  FH(i+6,i+2) = -J;
          FH(i+3,i+7) = -J;  FH(i+7,i+3) = -J;
          
        }

        FH(1,sizeH-3)=-J;  FH(sizeH-3,1)=-J;
        FH(2,sizeH-2)=-J;  FH(sizeH-2,2)=-J;  
        FH(3,sizeH-1)=-J;  FH(sizeH-1,3)=-J; 

      }
      else
        for(int i=0; i<sizeH-4; i=i+4){
          FH(i+1,i+5) = -2*J;
          FH(i+2,i+6) = -2*J;
          FH(i+3,i+7) = -2*J;

          FH(i+5,i+1) = -2*J;
          FH(i+6,i+2) = -2*J;
          FH(i+7,i+3) = -2*J;
        }

      return FH;
      }

      ComplexMatrix DrivenHamiltonian(){
        int sizeH=4*n_sites;  // # de estados y por tanto dimensión de la matriz

        RealMatrix DH(sizeH,sizeH,0.0);  

        // Términos diagonales (qZF)

        for(int j=0; j<sizeH ; j=j+4){
          DH(j+1,j+1)=q;
          DH(j+2,j+2)=q;
        }                          

        //cout << DH << endl;

      return DH;
      }
};

class Floquet : public Hamiltonian {
public:

    Floquet(int _n_sites, double _mu, double _J, double _q, double _T, int _harmonics)
    : Hamiltonian(_n_sites, _mu, _J, _q, _T), harmonics(_harmonics) {
    }

  ComplexMatrix DiagonalElements(int n){
    identity_matrix<double> id(4*n_sites);
    Identity = id;
    ComplexMatrix DE(sizeH,sizeH,0.0);


    auto Integrand=[this](double x) {
      return pow(cos(two_pi<double>()*x/T),2); 
    };

    Complex I=trapezoidal(Integrand,0.0,T)/T;

    DE=FreeHamiltonian()+I*DrivenHamiltonian()-n*two_pi<double>()*h*Identity/T;

    return DE;
  }

  ComplexMatrix NonDiagonalElements(int n, int m){
    ComplexMatrix NDE(sizeH,sizeH);

    auto Integrand=[this,m,n](double x) {
      return pow(cos(two_pi<double>()*x/T),2)*exp(-iu*x*(n-m)*two_pi<double>()/T); 
    };
    
    Complex I=trapezoidal(Integrand,0.0,T)/T;

    NDE=I*DrivenHamiltonian();

    return NDE;
  }

  ComplexMatrix FloquetHamiltonian(){

    ComplexMatrix FH(sizeH*(2*harmonics+1),sizeH*(2*harmonics+1),0.0);

    for(int i=-harmonics; i<=harmonics; i++){
      for(int j=-harmonics; j<=harmonics; j++){
        SubComplexMatrix FHij(FH,range((i+harmonics)*sizeH,(i+1+harmonics)*sizeH),range((j+harmonics)*sizeH,(j+1+harmonics)*sizeH));
        
        if(i+harmonics==j+harmonics){
          FHij=DiagonalElements(i+harmonics);
        }
        else{
          FHij=NonDiagonalElements(i+harmonics,j+harmonics);
        }
      }
    }
    return FH;
  }

public:
  int harmonics;
  int sizeH=4*n_sites;
  double h=1;
};

/*
class GraphquasienergiesMu{
  private:
    int n_sites;
    double J;
    double q;
    double alpha;
    double T;
    int harmonics;
    double mu0;
    double muf;
    double nPasos;
  public:
    graph(int _n_sites, double _J, double _q, double _alpha, double _T, int _harmonics,double _mu0, double _muf ,int _nPasos){
      n_sites=_n_sites;
      J=_J;
      q=_q;
      alpha=_alpha;
      T=_T;
      harmonics=_harmonics;
      mu0=_mu0;
      muf=_muf;
      nPasos=_nPasos;
    }
   void Eigenenergies(){
      ofstream kokito;
      std::ofstream outfile("EigenenergiesMu.dat", std::ofstream::trunc);
      kokito.open("EigenenergiesMu.dat");
      RealVector Evalues[int(nPasos)+1];
      RealVector EigenE;
      double sizeH=3*n_sites+1;
      for(int i=0;i<=nPasos;i++){
          TimeEvolution Hamil(n_sites,mu0+i*(muf-mu0)/nPasos,J,q,alpha,T,harmonics);
          Evalues[i]=Hamil.quasienergies();
      }
      for(int i=0;i< sizeH*(2.*harmonics+1.);i++){
        for(int j=0;j<=nPasos;j++){
          kokito<<Evalues[j](i)<<" "<<mu0+j*(muf-mu0)/nPasos<<endl;
          }
          kokito<<endl;
          kokito<<endl;
        }
      kokito.close();
    }
};

class GraphquasienergiesJ{
  private:
    int n_sites;
    double mu;
    double q;
    double alpha;
    double T;
    int harmonics;
    double J0;
    double Jf;
    double nPasos;
  public:
    graph2(int _n_sites, double _mu, double _q, double _alpha, double _T, int _harmonics,double _J0, double _Jf ,int _nPasos){
      n_sites=_n_sites;
      mu=_mu;
      q=_q;
      alpha=_alpha;
      T=_T;
      harmonics=_harmonics;
      J0=_J0;
      Jf=_Jf;
      nPasos=_nPasos;
    }
   void Eigenenergies(){
      ofstream kokito;
      std::ofstream outfile("EigenenergiesJ.dat", std::ofstream::trunc);
      kokito.open("EigenenergiesJ.dat");
      RealVector Evalues[int(nPasos)+1];
      RealVector EigenE;
      double sizeH=3*n_sites+1;
      for(int i=0;i<=nPasos;i++){
          TimeEvolution Hamil(n_sites,mu,J0+i*(Jf-J0)/nPasos,q,alpha,T,harmonics);
          Evalues[i]=Hamil.quasienergies();
      }
      for(int i=0;i< sizeH*(2.*harmonics+1.);i++){
        for(int j=0;j<=nPasos;j++){
          kokito<<Evalues[j](i)<<" "<<J0+j*(Jf-J0)/nPasos<<endl;
          }
          kokito<<endl;
          kokito<<endl;
        }
      kokito.close();
    }
};
*/

class TimeEvolution: public Floquet{
public:

    TimeEvolution(int _n_sites, double _mu, double _J, double _q, double _T, int _harmonics)
    : Floquet(_n_sites, _mu, _J, _q, _T, _harmonics){
      ComplexMatrix H(sizeH*(2*harmonics+1),sizeH*(2*harmonics+1));
      V.resize(sizeH*(2*harmonics+1),sizeH*(2*harmonics+1));
      E.resize(sizeH*(2*harmonics+1));

      H=FloquetHamiltonian();


      diag(H,E,V);
    }

    ComplexMatrix Vectors(){
      return V;
    }

    RealVector quasienergies(){
      return E;
    }

    ComplexVector Coefficients(ComplexVector v0){
      ComplexVector Coef(sizeH);
      ComplexMatrix V=Vectors();

      for(int j=0; j<sizeH; j++){

        ComplexVector Counter0(sizeH,0.0);
        ComplexVector v=column(V,j+sizeH*harmonics);

        for(int i=0; i<2*harmonics+1; i++){

          SubComplexVector fc(v, range(i*sizeH,(i+1)*sizeH));

          Counter0=Counter0+fc;
        }
        Coef(j)=inner_prod(conj(Counter0),v0);

      }

      double sum=0;

      for(int i=0; i<sizeH; i++){
        sum=sum+pow(abs(Coef(i)),2);
      }
      cout<<"Los coeficientes suman: "<< sum <<endl;
      
      return Coef;
    }

    ComplexVector TimeEvolutionState(ComplexVector Coef, double t){
      ComplexVector FinalState(sizeH);
      ComplexVector FloquetStates[sizeH];
      ComplexMatrix V=Vectors();
      RealVector E=quasienergies();
      
      for(int j=0; j<sizeH; j++){
        ComplexVector Counter(sizeH,0.0);
        ComplexVector v=column(V,j+sizeH*harmonics);

        for(int i=0;i<2*harmonics+1;i++){

          SubComplexVector fc(v, range(i*sizeH,(i+1)*sizeH));

          Counter=Counter+fc*exp(-iu*(i-harmonics)*two_pi<double>()*t/T);
       
        }

        FloquetStates[j]=Coef(j)*exp(-iu*E(j+sizeH*harmonics)*t/h)*Counter;
      
      }


      for(int i=0;i<sizeH;i++){
        FinalState=FinalState+FloquetStates[i];
      }

      return FinalState;
    }
  private:
    int sizeH=4*n_sites;
    double h=1;
    ComplexMatrix V;
    RealVector E;
  };

//---------------------------------------------------------------------------------------------------------------------------

int main () 
{
  //----------------------------------------- System Variables ------------------------------------------------------------

  int n_sites=30;   // # Sitios
  int harmonics=1;  // # de Armónicos
  double mu=1;    // Chemical Potential
  double J=1;     // Hopping
  double q=1;     // qZF
  double T=0.001;     // Periodo 
  Complex iu(0,1); // Unidad Imaginaria



  ComplexVector state0(4*n_sites,0.0); // Estado inicial
  
  /*state0(0)=1/sqrt(2);
  
  for(int k=1; k<n_sites+1 ; k++){
    state0(3*k-2)=1/sqrt(n_sites)*exp(-iu*k*two_pi<double>()/n_sites);
  } 

  for(int k=1; k<n_sites-1 ; k++){
    state0(4*k)=0;
  } 
  */

  //state0(5)=1/sqrt(4);
  state0(1)=1;
  
  cout<< "\n" << endl; 
  cout<< "Estado inicial: " << state0 << "\n"<< endl;


  double tf=10;
  double dt=0.0001;
  double npasos=tf/dt;
  

  TimeEvolution OPBH(n_sites,mu,J,q,T,harmonics);

  ComplexVector Coef(4*n_sites);
  Coef=OPBH.Coefficients(state0);

  cout << "\n" << endl;
  //------------------------------------------Time Evolution for one state-----------------------------------------------------------
  
  /*
  ComplexVector statef(4*n_sites,0.0); // Estado a proyectar

  statef(1)=1;

  cout<< "Estado a proyectar: " << statef << endl;
  
  ofstream kokito;
  kokito.open("Prob.dat"); 
  //kokito << fixed << setprecision(15);

  double Prob[int(npasos)];

  for(int i=0; i<=int(npasos); i++){

    ComplexVector v=OPBH.TimeEvolutionState(Coef,dt*i);

    Prob[i]=pow(abs(inner_prod(conj(statef),v)),2);
 
    kokito<<Prob[i]<<" "<<dt*i<<endl;
  }

  kokito.close();
  */

  //-------------------------------------Observables--------------------------------------------------------
  
  
  ComplexMatrix Sz(4*n_sites,4*n_sites,0.0);
  ComplexMatrix Sx(4*n_sites,4*n_sites,0.0);
  ComplexMatrix Sy(4*n_sites,4*n_sites,0.0);
  ComplexMatrix Q1(4*n_sites,4*n_sites,0.0);
  ComplexMatrix Q2(4*n_sites,4*n_sites,0.0);
  ComplexMatrix Q3(4*n_sites,4*n_sites,0.0); // Observable
  ComplexMatrix Q4(4*n_sites,4*n_sites,0.0);
  ComplexMatrix Q5(4*n_sites,4*n_sites,0.0);
  ComplexMatrix n0(4*n_sites,4*n_sites,0.0);
  ComplexMatrix n1(4*n_sites,4*n_sites,0.0);
  ComplexMatrix n2(4*n_sites,4*n_sites,0.0);

  //Densidad Magnetizacion
  for(int i=0; i<4*n_sites; i=i+4){
    Sz(i+1,i+1)=1;

    Sz(i+2,i+2)=-1;
  }
  //Sz /= n_sites;
  

  //Magnetizacion Transversal

  for(int i=0; i<4*n_sites; i=i+4){
    Sx(i+3,i+1)=1/sqrt(2);
    Sx(i+1,i+3)=1/sqrt(2);
    
    Sx(i+3,i+2)=1/sqrt(2);
    Sx(i+2,i+3)=1/sqrt(2);

  }
  //Magnetizacion Transversal

  for(int i=0; i<4*n_sites; i=i+4){
    Sy(i+3,i+1)=-iu/sqrt(2);
    Sy(i+1,i+3)=iu/sqrt(2);

    
    Sy(i+3,i+2)=iu/sqrt(2);
    Sy(i+2,i+3)=-iu/sqrt(2);

  }

  //Q1
  for(int i=0; i<4*n_sites; i=i+4){
    Q1(i+3,i+1)=-iu/sqrt(2);
    Q1(i+1,i+3)=iu/sqrt(2);
    
    Q1(i+3,i+2)=-iu/sqrt(2);
    Q1(i+2,i+3)=iu/sqrt(2);

  }

  for(int i=0; i<4*n_sites; i=i+4){
    Q2(i+3,i+1)=1/sqrt(2);
    Q2(i+1,i+3)=1/sqrt(2);
   
    Q2(i+3,i+2)=-1/sqrt(2);
    Q2(i+2,i+3)=-1/sqrt(2);

  }

  for(int i=0; i<4*n_sites; i=i+4){
    Q3(i+2,i+1)=-iu/sqrt(2);
    Q3(i+1,i+2)=iu/sqrt(2);

  }

  for(int i=0; i<4*n_sites; i=i+4){
    Q4(i+2,i+1)=1;
    Q4(i+1,i+2)=1;

  }
  //Densidad Cuadrupolar

  for(int i=0; i<4*n_sites; i=i+4){
    Q5(i+1,i+1)=1/sqrt(3);

    Q5(i+2,i+2)=1/sqrt(3);

    Q5(i+3,i+3)=-2/sqrt(3);
  }
  //q /= n_sites;

  //Numero de partículas en 0
  for(int i=0; i<4*n_sites; i=i+4){
    n0(i+3,i+3)=1;
  }

  //Numero de partículas en 1
  for(int i=0; i<4*n_sites; i=i+4){
    n1(i+1,i+1)=1;
  }

  //Numero de partículas en -1
  for(int i=0; i<4*n_sites; i=i+4){
    n2(i+2,i+2)=1;
  }

  

  //cout<< "Observable: " << q << "\n"<< endl;

  double ExpectedValue[int(npasos)];

  ofstream kokito2;
  ofstream kokito3;
  ofstream kokito4;
  ofstream kokito5;
  ofstream kokito6;
  ofstream kokito7;
  ofstream kokito8;
  ofstream kokito9;
  ofstream kokito10;
  ofstream kokito11;
  ofstream kokito12;

  kokito2.open("Densidad_quadrupolar.dat");
  kokito3.open("Densidad_Magnetica.dat");
  kokito4.open("N_Particles_0.dat");
  kokito5.open("N_Particles_1.dat");
  kokito6.open("N_Particles_2.dat");
  kokito7.open("N_Particles_Total.dat");
  kokito8.open("Magnetization_Transversal_x.dat");
  kokito9.open("Magnetization_Transversal_y.dat");
  kokito10.open("Suma_Magnetica.dat");
  kokito11.open("Suma_Quadrupolar.dat");
  kokito12.open("Suma_Total.dat");

  kokito2 << fixed << setprecision(10);
  kokito3 << fixed << setprecision(10);
  kokito4 << fixed << setprecision(10);
  kokito5 << fixed << setprecision(10);
  kokito6 << fixed << setprecision(10);
  kokito7 << fixed << setprecision(10);
  kokito8 << fixed << setprecision(10);
  kokito9 << fixed << setprecision(10);
  kokito10 << fixed << setprecision(10);
  kokito11 << fixed << setprecision(10);
  kokito12 << fixed << setprecision(10);

  for(int i=0; i<=int(npasos); i++){

    ComplexVector v=OPBH.TimeEvolutionState(Coef,dt*i);

    ExpectedValue[i]=real(inner_prod(conj(v),prod(Q5,v)));

    //cout << inner_prod(conj(v),prod(q,v)) << endl;

    kokito2<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=real(inner_prod(conj(v),prod(Sz,v)));

    kokito3<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=real(inner_prod(conj(v),prod(n0,v)));

    kokito4<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=real(inner_prod(conj(v),prod(n1,v)));

    kokito5<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=real(inner_prod(conj(v),prod(n2,v)));

    kokito6<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=real(inner_prod(conj(v),prod(n0,v)))+real(inner_prod(conj(v),prod(n1,v)))+real(inner_prod(conj(v),prod(n2,v)));

    kokito7<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=pow(abs(inner_prod(conj(v),prod(Sx,v))),2);

    kokito8<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=pow(abs(inner_prod(conj(v),prod(Sy,v))),2);

    kokito9<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=pow(real(inner_prod(conj(v),prod(Sx,v))),2)+pow(real(inner_prod(conj(v),prod(Sy,v))),2)+pow(real(inner_prod(conj(v),prod(Sz,v))),2);

    kokito10<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=pow(real(inner_prod(conj(v),prod(Q1,v))),2)+pow(real(inner_prod(conj(v),prod(Q2,v))),2)+pow(real(inner_prod(conj(v),prod(Q3,v))),2)+pow(real(inner_prod(conj(v),prod(Q4,v))),2)+pow(real(inner_prod(conj(v),prod(Q5,v))),2);

    kokito11<<ExpectedValue[i]<<" "<<dt*i<<endl;

    ExpectedValue[i]=pow(real(inner_prod(conj(v),prod(Sx,v))),2)+pow(real(inner_prod(conj(v),prod(Sy,v))),2)+pow(real(inner_prod(conj(v),prod(Sz,v))),2)+pow(real(inner_prod(conj(v),prod(Q1,v))),2)+pow(real(inner_prod(conj(v),prod(Q2,v))),2)+pow(real(inner_prod(conj(v),prod(Q3,v))),2)+pow(real(inner_prod(conj(v),prod(Q4,v))),2)+pow(real(inner_prod(conj(v),prod(Q5,v))),2);

    kokito12<<ExpectedValue[i]<<" "<<dt*i<<endl;
  }

  kokito2.close();
  kokito3.close();
  kokito4.close();
  kokito5.close();
  kokito6.close(); 
  kokito7.close(); 
  kokito8.close();
  kokito9.close();
  kokito10.close(); 
  kokito11.close();
  kokito12.close();
   
  
  //-------------------------------------Espectral Analysis (quasienergies)-----------------------------------

  /*
  
  Floquet OPBH(n_sites,mu,J,q,alpha,T,harmonics);

  //cout<<OPBH.FreeHamiltonian()<<endl;
  //cout<<OPBH.DrivenHamiltonian()<<endl;
  cout<<OPBH.Energies()<<"\n"<<endl;

  //(sitios, J, q, T, #armonicos, mu_0,mu_f, #pasos)

  
  double mu0=-1.; // Chemical potential inicial
  double muf=1.; // Chemical potential final 
  double nPasos=300; // # de pasos
  
  GraphquasienergiesMu OPBH_mu(n_sites,J,q,alpha,T,harmonics,mu0,muf,nPasos);

  OPBH_mu.Eigenenergies();

  cout<< "quasienergias Calculadas-Mu"<< "\n" << endl;

  double J0=-1.; // Chemical potential inicia         l
  double Jf=1.; // Chemical potential final  

  GraphquasienergiesJ OPBH_J(n_sites,mu,q,alpha,T,harmonics,J0,Jf,nPasos);

  OPBH_J.Eigenenergies();

  cout<< "quasienergias Calculadas-J"<<endl;
 
  */

  return 0;
}