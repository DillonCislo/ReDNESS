/* ==============================================================================================
 *
 *  evaluate_elastic_energy.cpp
 *
 *  Evaluates the elastic energy of a non-Euclidean shell given a user-supplied
 *  real-space configuration and target geometry. Depends on CGAL.
 *
 *  by Dillon Cislo
 *  10/4/2019
 *
 *  This is a MEX-File for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <stdexcept>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../../ElasticEnergy/StretchOperator.h"
#include "../../ElasticEnergy/BendOperator.h"
#include "../../ElasticEnergy/HingeOperator.h"
#include "../../ElasticEnergy/FixedPointOperator.h"

#include "../../ElasticMesh/ElasticUpdater.h"

typedef CGAL::Simple_cartesian<double>              Kernel;
typedef Kernel::Point_3                             Point;
typedef Kernel::Vector_3                            Vector;

typedef CGAL::Polyhedron_3<Kernel, ElasticItems>    Polyhedron;
typedef Polyhedron::Vertex_handle                   Vertex_handle;
typedef Polyhedron::Vertex_iterator                 Vertex_iterator;
typedef Polyhedron::HalfedgeDS                      HalfedgeDS;

typedef Eigen::VectorXd                 VectorXd;
typedef Eigen::MatrixXd                 MatrixXd;
typedef Eigen::VectorXi                 VectorXi;
typedef Eigen::MatrixXi                 MatrixXi;

namespace PMP = CGAL::Polygon_mesh_processing;


///
/// Main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

  // ------------------------------------------------------------------------------------
  // INPUT PROCESSING
  // ------------------------------------------------------------------------------------
  
  // Check for proper number of arguments
  if ( nrhs != 9 ) {
    mexErrMsgIdAndTxt( "MATLAB:evaluate_elastic_energy:nargin",
        "EVALUATE_ELASTIC_ENERGY requires nine input arguments." );
  } else if ( nlhs != 3 ) {
    mexErrMsgIdAndTxt( "MATLAB:evaluate_elastic_energy:nargou",
        "EVALUATE_ELASTIC_ENERGY requires three output arguments." );
  }

  double *face = mxGetPr( prhs[0] );          // The face connectivity list
  std::size_t numFaces = mxGetM( prhs[0] );   // The number of faces
  std::size_t sizeFaces = mxGetN( prhs[0] );  // The number of vertices in a single face

  // Check that the mesh is a triangulation
  if ( sizeFaces != 3 ) {
    mexErrMsgIdAndTxt( "MATLAB:evaluate_elastic_energy:face_degree",
        "Faces must be elements of a triangulation." );
  }

  double *vertex = mxGetPr( prhs[1] );        // The vertex coordinate list
  std::size_t numVertex = mxGetM( prhs[1] );  // The number of vertices
  std::size_t dim = mxGetN( prhs[1] );        // The dimension of the vertex coordinates

  // Check the dimensions of the vertex coordinates
  if ( dim != 3 ) {
    mexErrMsgIdAndTxt( "MATLAB:evaluate_elastic_energy:vertex_dim",
        "Vertex coordinates must be 3D." );
  }

  double *m_targetLengths = mxGetPr( prhs[2] ); // The list of target edge lengths
  double *m_targetAngles = mxGetPr( prhs[3] );  // The list of target bend angles

  int *m_v1_ID = (int*) mxGetData( prhs[4] );     // The start vertex ID of each edge
  int *m_v2_ID = (int*) mxGetData( prhs[5] );     // The end vertex ID of eage edge
  std::size_t numEdges = mxGetM( prhs[5] ); // The number of edges in the triangulation

  double h = *mxGetPr( prhs[6] );   // The thickness of the shell
  double nu = *mxGetPr( prhs[7] );  // Poisson's ratio
  double Y = *mxGetPr( prhs[8] );   // Young's modulus


  // Re-format input parameters as Eigen-style vectors for processing -------------------
  
  VectorXi v1_ID = Eigen::Map<VectorXi>( m_v1_ID, numEdges );
  VectorXi v2_ID = Eigen::Map<VectorXi>( m_v2_ID, numEdges );

  VectorXd targetLengths = Eigen::Map<VectorXd>( m_targetLengths, numEdges );
  VectorXd targetAngles = Eigen::Map<VectorXd>( m_targetAngles, numEdges );

  // Create and populate the polyhedral mesh ---------------------------------------------
  
  // Create vector of 3D point objects
  std::vector<Point> points;
  points.reserve( numVertex );
  for( int i = 0; i < numVertex; i++ ) {

    points.push_back( Point( vertex[i],
          vertex[i+numVertex],
          vertex[i+(2*numVertex)] ) );

  }

  // Create vector of polygon objects
  std::vector< std::vector<std::size_t> > polygons;
  polygons.reserve( numFaces );
  for( int i = 0; i < numFaces; i++ ) {

    std::vector<std::size_t> currentPolygon;
    currentPolygon.reserve( sizeFaces );
    for( int j = 0; j < sizeFaces; j++ ) {

      // NOTE: We subtract 1 from the index to account for
      // MATLAB's 1-indexed array structures
      double index = face[i+(j*numFaces)]-1.0;
      currentPolygon.push_back( (std::size_t) index );

    }

    polygons.push_back( currentPolygon );

  }

  // Populate the mesh
  Polyhedron P;
  PMP::orient_polygon_soup( points, polygons );
  PMP::polygon_soup_to_polygon_mesh( points, polygons, P );

  // Set mesh geometry ---------------------------------------------------------------------
  
  ElasticUpdater EU;

  EU.assignVertexID( P );
  EU.assignMajorEdges( P );

  EU.updateTargetGeometry( P, v1_ID, v2_ID, targetLengths, targetAngles );

  VectorXd x = Eigen::Map<VectorXd>( vertex, dim*numVertex );
  EU.updateCurrentGeometry( P, x );

  // ----------------------------------------------------------------------------------------
  // CALCULATE ELASTIC ENERGY
  // ----------------------------------------------------------------------------------------
  
  StretchOperator SO( h, nu, Y );
  BendOperator BO( h, nu, Y );

  double Estretch = SO( P );
  double Ebend = BO( P );

  double Etotal = Estretch + Ebend;

  // ---------------------------------------------------------------------------------------
  // OUTPUT PROCESSING
  // ---------------------------------------------------------------------------------------
  
  plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
  double *ETO = mxGetPr( plhs[0] );
  ETO[0] = Etotal;

  plhs[1] = mxCreateDoubleMatrix( 1, 1, mxREAL );
  double *ESO = mxGetPr( plhs[1] );
  ESO[0] = Estretch;

  plhs[2] = mxCreateDoubleMatrix( 1, 1, mxREAL );
  double *EBO = mxGetPr( plhs[2] );
  EBO[0] = Ebend;

};
