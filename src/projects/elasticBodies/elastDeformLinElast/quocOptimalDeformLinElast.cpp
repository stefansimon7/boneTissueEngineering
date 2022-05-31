#include <quocConfigurators.h>
#include <quocDefines.h>

#include "quocOptimalDeformSolverLinElast.h"

using namespace quocFE;

#define _USEQUOC2D
// #define _USEQUOC3D

#ifdef _USEQUOC2D
static const int dimDomain=2;
typedef Quoc2DDataTypeContainer                                                                                       DataTypeContainer;
typedef QuocConfigurator2D<DataTypeContainer>                                                                         ConfiguratorType;
#endif

#ifdef _USEQUOC3D
static const int dimDomain=3;
typedef Quoc3DDataTypeContainer                                                                                        DataTypeContainer;
typedef QuocConfigurator3D<DataTypeContainer>                                                                          ConfiguratorType;
#endif

typedef typename DataTypeContainer::RealType                                                                            RealType;
typedef typename DataTypeContainer::VectorType                                                                          VectorType;
typedef QuocMaterialOptimizationConfigurator< ConfiguratorType >                                                        MatOptConfigurator;
typedef typename DataTypeContainer::ParameterParserType                                                                 ParameterParserType;
typedef typename DataTypeContainer::PointType                                                                           PointType;
typedef typename ConfiguratorType::SparseMatrixType                                                                     SparseMatrixType;
typedef typename ConfiguratorType::MaskType                                                                             MaskType;
typedef typename ConfiguratorType::InitType                                                                             MeshType;


int main () {
    aol::consoleOutputStartProgramm( "Find Optimal Deformation (linear elasticity)" );
    
    //! start watch
    auto startTime = std::chrono::high_resolution_clock::now();
    
    //! Parser
    ParameterParserType parser( "../../../../parser/shapeDesignBones/ElastShapeOpt.ini", "../../../../parser/counter.txt", "/LinElastOptDeformation" );
    
    //! save to log-file
    aol::AdditionalOutputToFile pAddOut ( aol::strprintf ( "%s/log.txt", parser.template get<string> ( "saving.saveDirectory" ).c_str () ) );
   
    //! Initialize mesh and configurator 
#ifdef _USEQUOC2D
    MeshType mesh ( parser.template get<int>("InputMesh.Nx"),parser.template get<int>("InputMesh.Ny"), parser.template get<int>("InputMesh.lx"),parser.template get<int>("InputMesh.ly") );
#endif
    
#ifdef _USEQUOC3D
    MeshType mesh ( parser.template get<int>("InputMesh.Nx"),parser.template get<int>("InputMesh.Ny"),parser.template get<int>("InputMesh.Nz"),
                    parser.template get<int>("InputMesh.lx"),parser.template get<int>("InputMesh.ly"),parser.template get<int>("InputMesh.lz") );
#endif
    ConfiguratorType conf ( mesh );
    MatOptConfigurator matOpConf ( parser, conf );
    
    cout << endl << "We have " << endl 
         << mesh.getNumElements()  << " Elements" << endl
         << mesh.getNumVertices() << " NumVertices" << endl
         << "dim domain = " << dimDomain << endl;
    
    //initialize material
    VectorType material ( conf.getNumGlobalDofs() );
     for( int nodeIdx = 0; nodeIdx < mesh.getNumVertices(); ++nodeIdx ){
         material[nodeIdx] = 0.0;
    }
      
    QuocHandler<ConfiguratorType> quocHandler ( parser, conf );
    
    //FORCE
    VectorType rhs_force ( dimDomain * conf.getNumGlobalDofs() );
    quocCreateForce<ConfiguratorType> ( parser, conf, rhs_force, quocHandler.getDirichletMask() );
    
    //! solve
    quocOptimalDeformSolverLinElast<MatOptConfigurator> linElastProblem( parser, matOpConf, material, quocHandler, rhs_force );
    
    
    //PLOT
    const VectorType &solution  = linElastProblem.getSolutionDisplacement();
    MeshType meshToPlot ( mesh );
    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        PointType coords = mesh.getVertex( i );
        for( int comp = 0; comp < dimDomain; ++comp ) coords[comp] += solution[i + comp * conf.getNumGlobalDofs()];
        meshToPlot.setVertex( i, coords );
    }
    MeshWithData<MeshType> meshSaver ( meshToPlot );
    meshSaver.saveAsVTKPOLYDATA( aol::strprintf ( "%s/solutionDisp.vtk", parser.template get<std::string> ( "saving.saveDirectory" ).c_str () ) );
    quocHandler.plotWithMaterial ( linElastProblem.getSolutionDisplacement(), material, "solDisp" );
  
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << "duration = " << diff.count() << " sec" << endl;

  return ( EXIT_SUCCESS );
}
