#ifndef PARAMS_CLASS_H
#define PARAMS_CLASS_H




class ParamsClass
{
public:
  double beta;
  double dipfac; 
  double qO;
  double qH;

  double cOO;
  double cOH;
  
  double bOO;
  double bOH;
  double alpha_O;
  double raggio;


  double D1_OO;
  double D1_OH;
  double D1_HH;

  double D2_OO;
  double D2_OH;
  double D2_HH;

  double gamma1_OO;
  double gamma1_OH;
  double gamma1_HH;

  double gamma2_OO;
  double gamma2_OH;
  double gamma2_HH;


  double r1_OO;
  double r1_OH;
  double r1_HH;



  double r2_OO;
  double r2_OH;
  double r2_HH;



  dVecp Box;
  dVecp BoxInv;
  void Init()
  {
    Box.vec[0]= 12.431071665375*angs2bohr; 
    Box.vec[1]= 12.431071665375*angs2bohr; 
    Box.vec[2]= 12.431071665375*angs2bohr; 
    beta=0.75;
    for (int dim=0;dim<NDIM;dim++)
      BoxInv.vec[dim]=1.0/Box.vec[dim];
      
    raggio=1.5; // 2.9655460304201733; 

    qO=-1.1995807E+00;
      //    qO=-1.1995809;
    qH=0.59979035;

    cOO=6.8676628;
    //    cOO=6.86766;
    //    cOH=-2.04454;
    cOH=-2.0445413;

    //    bOO=2.24850;
    bOO=2.2485039;
    
    //    bOH=3.92473;
    bOH=3.9247332;

    //    alpha_O=4.08675;
    alpha_O=4.0867573;
    //    dipfac=1.01977602741587826E-002;
    dipfac = 2.0/(3.0*pow(raggio,3.0)*sqrt(2.0*M_PI));


    D1_OO=-1.26367e-4;
    D1_OH=3.77059e-5;
    D1_HH=8.71664e-1;

    D2_OO=3.08029e-4;
    D2_OH=8.67779e-6;
    D2_HH=1.93115e-7;

    

    gamma1_OO=13.09317;
    gamma1_OH=15.20544;
    gamma1_HH=13.13611;

    gamma2_OO=13.96521;
    gamma2_OH=12.38136;
    gamma2_HH=16.13997;

    r1_OO=7.473237;
    r1_OH=3.08451;
    r1_HH=0.38486;

    r2_OO=7.03818;
    r2_OH=5.63316;
    r2_HH=8.12431;



  }

};


#endif
