/*-------------------------------------------------------------*/
/*  how to use the NOMAD library with a FORTRAN user function  */
/*  JCM: Original version basic_lib.cpp                        */
/*-------------------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

#define USE_SURROGATE false 

extern"C" {
  void bb_kleed_ ( double x[42] , double fx[1] );
}

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
public:
  My_Evaluator  ( const Parameters & p ) :
    Evaluator ( p ) {}

  ~My_Evaluator ( void ) {}

  bool eval_x ( Eval_Point   & x          ,
		const Double & h_max      ,
		bool         & count_eval   ) const {

    double xx[42];
    double fx[1];
    int    i;
    //    cout << "eval_x: values of x\n";

    for ( i = 0 ; i < 42 ; ++i ) {
      xx[i] = x[i].value();
      //      cout << i << "    " << xx[i] << "\n" ;
    }
    // call the FORTRAN routine:
    // cout << "Calling bb_kleed\n";
    bb_kleed_ ( xx , fx );

    for ( i = 0 ; i < 1 ; ++i )
      x.set_bb_output ( i , fx[i] );

    count_eval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
  }
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // display:
  Display out ( std::cout );
  out.precision ( DISPLAY_PRECISION_STD );

  // NOMAD initializations:
  begin ( argc , argv );

  try {

    // parameters creation:
    Parameters p ( out );
    int i;

    if ( USE_SURROGATE )
      p.set_HAS_SGTE ( true );

    p.set_DIMENSION (42);             // number of variables
    vector<bb_output_type> bbot (1); // definition of
    bbot[0] = OBJ;                 // output types
    p.set_BB_OUTPUT_TYPE ( bbot );

    Point x0(42);
    //    double xstart[42] = {-1.8757, -1.7941, -1.8067, -0.3861, 0.2472, -0.0461, 0.0690, 0.1874, 1.7112, 1.7350 
    //	      1.7378, 1.7467, 1.7751, 1.7897, 0.0, 3.1141, 3.0047, 6.2250, -4.0621, 1.2552, 3.6738, 
    //	      -4.2907, 5.0398, 0.0, 5.0355, 2.4703, 2.5445, 2.4371, 0.0, 0.0, -3.0047, 1.2913, 6.2250, 1.2552,
    //			 1.2125, 3.7093, 0.0, 0.0, 5.0355, 5.0402, 0.0, 2.4371};

    double xstart[42] = {-1.87570000000000, -1.79410005000000, -1.80669999000000, -0.38609999000000, 
			 -0.25279999000000, -0.04610000000000, 0.06900000000000, 0.18740000000000, 
			 1.71120000000000, 1.73500001000000, 1.73780000000000, 1.74670005000000, 
			 1.77509999000000, 1.78970003000000, 
			 0.0, 3.11409998000000, 3.00469995000000, 6.22499990000000, 
			 3.93790007000000, 1.25520003000000, 3.67379999000000, 3.70930004000000, 
			 5.03980017000000, 0.0, 5.03550005000000, 2.47029996000000, 
			 2.54450011000000, 2.43709993000000, 
			 0.0, 0.0, -3.00469995000000, 1.29130006000000, 
			 6.22499990000000, 1.25520003000000, 1.21249998000000, 3.70930004000000, 
			 0.0, 0.0, 5.03550005000000, 5.04020023000000, 
			 0.0, 2.43709993000000};

    for ( i = 0 ; i < 42 ; ++i ) {
      x0[i] = xstart[i];  // cout << i << "    " << x0[i] << "\n" ;
    }

    p.set_X0 (x0);  // starting point

    Point lb (42), ub(42);
    for ( i = 0 ; i < 14 ; ++i ) {
      lb[i] = -3.0; lb[i+14] = -4.0; lb[i+28] = -4.0; 
      ub[i] =  3.0; ub[i+14] = 10.0; ub[i+28] =  10.0; 
    }
    p.set_LOWER_BOUND (lb);
    p.set_UPPER_BOUND (ub);
    p.set_INITIAL_MESH_SIZE(1.0);

    p.set_MAX_BB_EVAL (2000);     // the algorithm terminates after
                                 // 500 black-box evaluations
    // p.set_TMP_DIR ("/tmp");      // repertory for
                                   // temporary files

    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev   ( p );

    // algorithm creation and execution:
    Mads mads ( p , &ev );
    mads.run();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  Slave::stop_slaves ( out );
  end();

  return EXIT_SUCCESS;
}
