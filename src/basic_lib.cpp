/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
// using namespace NOMAD; avoids putting NOMAD:: everywhere

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/

extern "C" {
double eval_kleed_fort_(double x[], double fitval[]);
}
// c stub for fortan subroutine
//
//  double fx = 0;
//  fitval = fx;
//  return(fitval);
//}
NOMAD::Double eval_kleed(NOMAD::Eval_Point &x) {
  double xparm[42];
  double fitval[3];
  for ( int i = 0 ; i < 42 ; i++ ) 
    {
      xparm[i] = x[i].value();
    }
  eval_kleed_fort_(xparm, fitval);
  return fitval[0];
}


class My_Evaluator : public NOMAD::Evaluator {
public:
  My_Evaluator  ( const NOMAD::Parameters & p ) :
    NOMAD::Evaluator ( p ) {}

  ~My_Evaluator ( void ) {}

  bool eval_x ( NOMAD::Eval_Point   & x          ,
		const NOMAD::Double & h_max      ,
		bool                & count_eval   ) const
	{

	  NOMAD::Double kx; 
	  int DIR = 0;
	  int RANK = 0;
	  int NMAX = 14;
	  int NDIM = 3;

	  kx = eval_kleed(x);
	  x.set_bb_output  ( 0 , kx  ); // objective value
	  x.set_bb_output  ( 1 , 0); // constraint 1

	  count_eval = true; // count a black-box evaluation
	  
	  return true;       // the evaluation succeeded
	}
	
  
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // display:
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  try {

    // NOMAD initializations:
    NOMAD::begin ( argc , argv );

    NOMAD::RNG::set_seed(12345);
	  

    // parameters creation:
    NOMAD::Parameters p ( out );

    p.set_DIMENSION (42);             // number of variables

    vector<NOMAD::bb_output_type> bbot (2); // definition of
    bbot[0] = NOMAD::OBJ;                   // output types
    bbot[1] = NOMAD::PB;
    p.set_BB_OUTPUT_TYPE ( bbot );

    //p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
    p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );

    p.set_X0 ( NOMAD::Point(42,0.0) );  // starting point

    p.set_LOWER_BOUND ( NOMAD::Point (42, -1.0 ) ); // all var. >= -1
    p.set_UPPER_BOUND ( NOMAD::Point (42,  1.0 ) ); // all var. <=  1

    p.set_MAX_BB_EVAL (100);     // the algorithm terminates after
                                 // 100 black-box evaluations
	  
    // p.set_TMP_DIR ("/tmp");   // directory for temporary files

    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev   ( p );

    // algorithm creation and execution:
    NOMAD::Mads mads ( p , &ev );
    mads.run();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  NOMAD::Slave::stop_slaves ( out );
  NOMAD::end();

  return EXIT_SUCCESS;
}
