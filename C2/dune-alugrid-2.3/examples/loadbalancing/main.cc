/******************************************************************************

  ALUGrid - a library providing a mesh manager supporting simplicial
  and hexahedral meshes, local grid adaptivity for use in parallel
  computations including dynamic load balancing.

  Copyright (C) 1998 - 2002 Bernhard Schupp
  Copyright (C) 1998 - 2002 Mario Ohlberger
  Copyright (C) 2004 - 2012 Robert Kloefkorn
  Copyright (C) 2005 - 2012 Andreas Dedner
  Copyright (C) 2010 - 2012 Martin Nolte

  The DUNE ALUGrid module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The ALUGrid library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

******************************************************************************/

/******************************************************************************

  ALUGrid - a library providing a mesh manager supporting simplicial
  and hexahedral meshes, local grid adaptivity for use in parallel
  computations including dynamic load balancing.

  Copyright (C) 1998 - 2002 Bernhard Schupp
  Copyright (C) 1998 - 2002 Mario Ohlberger
  Copyright (C) 2004 - 2012 Robert Kloefkorn
  Copyright (C) 2005 - 2012 Andreas Dedner
  Copyright (C) 2010 - 2012 Martin Nolte

  The DUNE ALUGrid module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The ALUGrid library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

******************************************************************************/

/** include config file generated by configure
 *  (i.e., know what grids are present, etc)
 *  this should always be included first */
#include <config.h>
/** standard headers **/
#include <iostream>
/** dune (mpi, field-vector and grid type for dgf) **/
#include <dune/common/mpihelper.hh>     
#include <dune/common/fvector.hh>        
#include <dune/common/timer.hh>        

/** numerical scheme **/
#include "../piecewisefunction.hh"

/** adaptation scheme **/
#include "loadbalance.hh"
#include "adaptation.hh"

template <class Grid>
struct AssignRank
{
  AssignRank(const int rank) : rank_(rank) {}
  Dune::FieldVector<double,2> 
    initial(const Dune::FieldVector<double,Grid::dimensionworld> &) const
  {
    return Dune::FieldVector<double,2>(rank_);
  }
  private:
  int rank_;
};

// method
// ------
void method ( int startLevel, int maxLevel, const char* outpath )
{
  typedef Dune::GridSelector::GridType Grid;
  /* Grid construction ... */
  std::string name = "../dgf/unitcube3d.dgf" ;
  // create grid pointer and release to free memory of GridPtr
  Grid* gridPtr = Dune::GridPtr<Grid>(name).release() ;

  Grid &grid = *gridPtr;

#if HAVE_ZOLTAN && USE_ZOLTANLB
  typedef ZoltanLoadBalanceHandle<Grid> LoadBalancer;
#else
  typedef SimpleLoadBalanceHandle<Grid> LoadBalancer;
#endif
  LoadBalancer ldb(grid);

  {
    typedef Dune::LoadBalanceHandleIF< LoadBalancer > DataHandleInterface;
    grid.loadBalance( (DataHandleInterface&)(ldb) );
    // grid.loadBalance();
  }

  const bool verboseRank = grid.comm().rank() == 0 ;

  std::string outPath( outpath );
  const bool writeOutput = ( outPath != "none" ) ;

  /* ... some global refinement steps */
  if( verboseRank ) 
    std::cout << "globalRefine: " << startLevel << std::endl;
  grid.globalRefine( startLevel );

  /* get view to leaf grid */
  typedef Grid::Partition< Dune::Interior_Partition >::LeafGridView GridView;
  GridView gridView = grid.leafView< Dune::Interior_Partition >();

  /* construct data vector for solution */
  typedef PiecewiseFunction< GridView, Dune::FieldVector< double, 2 > > DataType;
  DataType solution( gridView );
  solution.initialize( AssignRank<Grid>(grid.comm().rank()) );

  /* create VTK writer for data sequqnce */
  Dune::VTKSequenceWriter< GridView > vtkOut( gridView, "solution", outPath, ".", Dune::VTK::nonconforming );
  if( writeOutput ) 
  {
    VTKData< DataType >::addTo( solution, vtkOut );
    VTKData< DataType >::addPartitioningData( grid.comm().rank(), vtkOut );
  }

  /* create adaptation method */
  typedef LeafAdaptation< Grid, LoadBalancer > AdaptationType;
  AdaptationType adaptation( grid, ldb );

  if( writeOutput ) 
  {
    /* output the initial grid and the solution */
    vtkOut.write( 0.0 );
  }

  /* final time for simulation */
  const double endTime = 1.;
  /* interval for saving data */
  const double saveInterval = 0.02;
  /* first point where data is saved */
  double saveStep = saveInterval;

  /* now do the time stepping */
  unsigned int step = 0;
  double time = 0.0;
  while ( time < endTime ) 
  {
    double dt = saveInterval;

    /* augment time */
    time += dt;
    ++step;

    /* check if data should be written */
    if( time >= saveStep )
    {
      if( writeOutput ) 
      {
        /* visualize with VTK */
        vtkOut.write( time );
      }
      /* set saveStep for next save point */
      saveStep += saveInterval;
    }

    adaptation( solution );

  }           

  if( writeOutput ) 
  {
    /* output final result */
    vtkOut.write( time );
  }

  // delete grid 
  delete gridPtr ;
}
/***************************************************
 ** main program with parameters:                 **
 ** 1) number of problem to use (initial data...) **
 ** 2) number of global refinement steps          **
 ** 3) maximal level to use during refinement     **
 ***************************************************/
int main ( int argc , char **argv )
try
{
  /* initialize MPI, finalize is done automatically on exit */
  Dune::MPIHelper &mpi = Dune::MPIHelper::instance( argc, argv );
  
#if HAVE_ZOLTAN 
  float version;
  int rc = Zoltan_Initialize(argc, argv, &version);
  if (rc != ZOLTAN_OK){
    printf("sorry zoltan did not initialize successfully...\n");
    MPI_Finalize();
    exit(0);
  }
#endif

  if( argc < 1 )
  {
    /* display usage */
    if( mpi.rank() == 0 )
      std::cout << "Usage: " << argv[ 0 ] << " [startLevel] [maxLevel]" << std::endl;
    return 0;
  }

  /* get level to use for computationa */
  const int startLevel = (argc > 1 ? atoi( argv[ 1 ] ) : 0);
  const int maxLevel = (argc > 2 ? atoi( argv[ 2 ] ) : startLevel);

  const char* path = (argc > 3) ? argv[ 3 ] : "./";
  method( startLevel, maxLevel, path );

  /* done */
  return 0;
}
catch( const std::exception &e )
{
  std::cout << "STL ERROR: " << e.what() << std::endl;
  return 1;
}
catch( const Dune::Exception &e )
{
  std::cout << "DUNE ERROR: " << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cout << "Unknown ERROR" << std::endl;
  return 1;
}
