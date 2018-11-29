void quatroTQnl(double *dhdt, double *X, double *U){
// Finção que implementa EDO do sistema quatro tanques

// Definindo as constatntes ---------------------------------------------------------------------------------

static const double   A1 = 28.0,  //  [=] cm^2
                      A3 = 28.0,  //  [=] cm^2
                      A2 = 32.0,  //  [=] cm^2
                      A4 = 32.0;  //  [=] cm^2

static const double a1 = 0.071, //  [=] cm^2
                    a3 = 0.071, //  [=] cm^2
                    a2 = 0.057, //  [=] cm^2
                    a4 = 0.057; //  [=] cm^2

static const double kc = 0.50;  //  [=] V/cm

static const double g  = 981.0; //  [=] cm/s^2

static const double k1    = 3.14, //  [=] cm^3/V.s
                    k2    = 3.29, //  [=] cm^3/V.s
                    gama1 = 0.43, //  
                    gama2 = 0.34; //  

// Explicitando estados
double  h1 = X[0], //  [=] cm
  		  h2 = X[1], //  [=] cm
  		  h3 = X[2], //  [=] cm
  		  h4 = X[3]; //  [=] cm

// Explicitando entradas
double	v1 = U[0],    // [=] V
  		  v2 = U[1];    // [=] V

dhdt[0] = -(a1/A1)*sqrt(2.0*g*h1) + (a3/A1)*sqrt(2.0*g*h3) + gama1*k1*v1/A1;
dhdt[1] = -(a2/A2)*sqrt(2.0*g*h2) + (a4/A2)*sqrt(2.0*g*h4) + gama2*k2*v2/A2;
dhdt[2] = -(a3/A3)*sqrt(2.0*g*h3) + (1.0-gama2)*k2*v2/A3;
dhdt[3] = -(a4/A4)*sqrt(2.0*g*h4) + (1.0-gama1)*k1*v1/A4;

}