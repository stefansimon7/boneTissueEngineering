#include "quocMaterialOptimizationSolver.h"

#include <quocConfigurators.h>
#include <quocDefines.h>

using namespace quocFE;

#define _USEQUOC2D
// #define _USEQUOC3D

#ifdef _USEQUOC2D
static const int dimDomain=2;
typedef Quoc2DDataTypeContainer                                                                                         DataTypeContainer;
typedef typename DataTypeContainer::RealType                                                                            RealType;
typedef QuocConfigurator2D<DataTypeContainer>                                                                       ConfiguratorType;
#endif

#ifdef _USEQUOC3D
static const int dimDomain=3;
typedef Quoc3DDataTypeContainer                                                                                         DataTypeContainer;
typedef typename DataTypeContainer::RealType                                                                            RealType;
typedef QuocConfigurator3D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

typedef typename DataTypeContainer::VectorType                                                                          VectorType;
typedef QuocMaterialOptimizationConfigurator< ConfiguratorType >                                                        MatOptConfigurator;
typedef typename DataTypeContainer::ParameterParserType                                                                 ParameterParserType;
typedef typename DataTypeContainer::PointType                                                                           PointType;
typedef typename ConfiguratorType::SparseMatrixType                                                                     SparseMatrixType;
typedef typename ConfiguratorType::MaskType                                                                             MaskType;
typedef typename ConfiguratorType::InitType                                                                             MeshType;


//--------------------------------------------------------
int main(int argc, char ** argv) {
    
    aol::consoleOutputStartProgramm( "Phase field optimization for Linear Elasticity" );
    
    //! start watch
    auto startTime = std::chrono::high_resolution_clock::now(); 

    ParameterParserType parser( "../../../../parser/shapeDesignBones/ElastShapeOpt.ini", "../../../../parser/counter.txt", "/MaterialOptLinElastQuoc" );
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

    quocMaterialOptimizationSolver<MatOptConfigurator> matOptProblem ( parser );
    matOptProblem.computeOptimalDesign( mesh );
    
    //! print elapsed time
    std::chrono::duration<RealType> diff = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << "duration = " << diff.count() << " sec" << endl;

  return 0;
}
