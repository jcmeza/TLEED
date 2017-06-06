/*-------------------------------------------------------------*/
/*  how to use the NOMAD library with a FORTRAN user function  */
/*  JCM: Original version basic_lib.cpp                        */
/*-------------------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

//#define USE_SURROGATE true
#define USE_SURROGATE false

extern"C" {
  void bb_tleed_ ( double x[42] , double fx[1] );
}

//extern"C" {
//  void bb_kleed_ ( double x[42] , double fx[1] );
//}

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
    }
    // call the FORTRAN routine:
    if (x.get_eval_type() == NOMAD::TRUTH ) {
      bb_tleed_ ( xx , fx );
    }
    else { // Surrogate
      // cout << "Calling bb_kleed\n";
      //      bb_kleed_ ( xx , fx );
    }

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
    int maxbb = 500;

    if ( USE_SURROGATE )
      p.set_HAS_SGTE ( true );

    p.set_DIMENSION (42);             // number of variables
    vector<bb_output_type> bbot (1); // definition of
    bbot[0] = OBJ;                 // output types
    p.set_BB_OUTPUT_TYPE ( bbot );

    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
    p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );
    p.set_DISPLAY_STATS ( "bbe  obj" );

    Point x0(42);
    //    double xstart[42] = {-1.8, -1.7, -1.8, -0.38, 0.24, -0.04, 0.06, 0.18, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7,
    //			 0.0, 3.1, 3.0, 6.2, 4.0, 1.2, 3.6, 4.2, 5.0, 0.0, 5.0, 2.4, 2.5, 2.4,
    //			 0.0, 0.0, -3.0, 1.2, 6.2, 1.2, 1.2, 3.7, 0.0, 0.0, 5.0, 5.0, 0.0, 2.4};

    //    double xstart[42] = { -2.088331226, -1.791994337, -1.795301675, -0.3909828025, -0.2518234275, -0.04281936035, 0.06375097656, 0.076071875, 1.741961719, 1.841079112, 1.767951367, 1.702266456, 1.829802749, 1.773052691, -0.2846679688, 5.132288457, 2.494171386, 5.822167869, 3.90324736, 1.001858355, 3.889864443, 3.755320548, 4.603826049, 0.8571166992, 5.059914112, 2.455758334, 2.593709705, 2.571117874, 0.00520324707, 0.4510498047, -2.663391356, 1.224878551, 5.51267385, 2.154430987, 1.109716777, 4.471629142, 0.8609619141, -0.3156738281, 4.756691456, 5.186440464, 0.8152313232, 3.449551102};
  
    double xstart[42] = {
      -1.8757, -1.7941, -1.8067, -0.3861, -0.2528, -0.0461, 0.0690, 0.1874, 1.7112, 1.735, 1.7378, 1.7467, 1.7751, 1.7897, 
      0.0,  3.1141, 3.0047, 6.2250, -4.0621, 1.2552, 3.6738, -4.2907, 5.0398, 0., 5.0355, 2.4703, 2.5445, 2.4371, 
      0.0, 0.0, -3.0047, 1.2913, 6.2250, 1.2552, 1.2125, 3.7093, 0., 0.,5.0355, 5.0402, 0., 2.4371};

    //cout << "nomad_tleed: Initializing x0\n";

    for ( i = 0 ; i < 42 ; ++i ) {
      x0[i] = xstart[i];
    }

    p.set_X0 (x0);  // starting point

    Point lb (42), ub(42);
    for ( i = 0 ; i < 14 ; ++i ) {
      lb[i] = -2.0; lb[i+14] = -5.0; lb[i+28] =  -4.0; 
      ub[i] =  2.0; ub[i+14] =  7.0; ub[i+28] =   7.0; 
    }
    p.set_LOWER_BOUND (lb);
    p.set_UPPER_BOUND (ub);
    p.set_INITIAL_MESH_SIZE(.1);

    p.set_MAX_BB_EVAL (maxbb);     // the algorithm terminates after
                                 // 500 black-box evaluations
    // p.set_TMP_DIR ("/tmp");      // repertory for
                                   // temporary files
    p.set_DISPLAY_DEGREE(2);
    p.set_SOLUTION_FILE("sol.txt");

    cout << "nomad_tleed: Initial x0\n";
    //    cout << "nomad_tleed: MAX_BB_EVAL = ", maxbb, "\n";
    for ( i = 0 ; i < 14 ; ++i ) {
      cout << i+1 << "    " << x0[i] << "    " << x0[i+14] << "    " << x0[i+28] << "\n" ;
    }

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
