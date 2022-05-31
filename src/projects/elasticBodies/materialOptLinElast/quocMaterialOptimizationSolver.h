#ifndef __QUOCMATERIALOPTIMIZATIONSOLVER_H
#define __QUOCMATERIALOPTIMIZATIONSOLVER_H

#include <derivativeTest.h>
#include <ipoptBoxConstraintSolver.h>
#include <ipoptNonlinearConstraintSolver.h>
#include <LineSearchMethods.h>
#include <quocHandler.h>
#include <SolverInfo.h>

#include "../elastDeformLinElast/quocOptimalDeformSolverLinElast.h"
#include "quocMaterialOptimizationComplianceEnergies.h"
#include "quocMaterialOptimizationEnergies.h"

using namespace quocFE;

//#define DEBUGINTERPOLATIONOPERATORADAPTIVE

template<typename MatOptConfigurator>
class quocMaterialOptimizationSolver {
  
  typedef typename MatOptConfigurator::ConfiguratorType         ConfiguratorType;
  typedef typename ConfiguratorType::RealType                   RealType;
  typedef typename ConfiguratorType::MaskType                   MaskType;
  typedef typename ConfiguratorType::InitType                   MeshType;
  typedef typename ConfiguratorType::DTContainer                DataTypeContainer;
  typedef typename DataTypeContainer::PointType                 PointType;
  typedef typename ConfiguratorType::VectorType                 VectorType;
  typedef typename ConfiguratorType::SparseMatrixType         SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType       ParameterParserType;

  const ParameterParserType &_parser;
  const int _maxIterations;
  const RealType _tolerance;
  mutable RealType _lastVolumeEnergy, _lastAreaEnergy, _lastJphysEnergy; // _lastTotalEnergy;

public:
  quocMaterialOptimizationSolver ( const ParameterParserType &Parser ) :
  _parser( Parser ),
  _maxIterations ( Parser.template get<int> ("MaterialOptimization.maxIterations") ),
  _tolerance( Parser.template get<RealType> ("MaterialOptimization.tolerance") ) { } 

  
  void getLastEnergies( RealType &Compliance, RealType &Area, RealType &Interface ) const{ 
      Compliance = _lastJphysEnergy; Area = _lastVolumeEnergy; Interface = _lastAreaEnergy;
  }
  
  void computeOptimalDesign( const MeshType & mesh ) const {
                                 
    aol::consoleOutput( "Compute optimal design " );
    
    ConfiguratorType conf ( mesh );
    QuocHandler<ConfiguratorType> quocHandler ( _parser, conf );
    
    //! initialize force
    VectorType rhs_force ( conf.dimDomain * conf.getNumGlobalDofs() );
    quocCreateForce<ConfiguratorType> ( _parser, conf, rhs_force, quocHandler.getDirichletMask() );
    
    //! initialize material
    VectorType initMaterial ( conf.getNumGlobalDofs() );
    quocHandler.switchMaterialType ( initMaterial );
//     material.setZero();
    VectorType solutionMaterial ( initMaterial );
    
    //!initialize disp
    VectorType solDisp ( conf.dimDomain * conf.getNumGlobalDofs( ) ); solDisp.setZero();
    
    //!initialize eps_area
    RealType eps_area = mesh.getInterfaceWith() * _parser.template get<RealType>( "MaterialOptimization.eps_factor" );
    
    //initialize MaterialOptimizationConfigurator
    MatOptConfigurator matOptConf ( _parser, conf, eps_area );
    
    //output
    cout << endl << "We have " << endl 
         << mesh.getNumElements()  << " Elements" << endl
         << mesh.getNumVertices() << " NumVertices" << endl
         << eps_area << " interface width" << endl << endl;
    
    //! initialize constraint
    cout << "1" << endl;
    quocOptimalDeformSolverLinElast<MatOptConfigurator> OptDeformFinder ( _parser, matOptConf, initMaterial, quocHandler, rhs_force );
    
    cout << "2" << endl;
    VectorType initMaterialScaledTo01( initMaterial.size() );
    for( int i=0; i<initMaterial.size(); i++ ) initMaterialScaledTo01[i] = 0.5 * ( initMaterial[i] + 1.0 );
    quocHandler.plotWithMaterial ( OptDeformFinder.getSolutionDisplacement(), initMaterialScaledTo01, "solutionDeformationForInitMaterial" );
    
    //! initialize adjoint problem
    cout << "3" << endl;
    MaterialOptimizationAdjointProblem<MatOptConfigurator> adjointProblem ( OptDeformFinder );
    
    //! construct energy operators
    cout << "4" << endl;
    MaterialOptimizationEnergyOp<MatOptConfigurator>  totalEnergyOp ( matOptConf, OptDeformFinder, adjointProblem );
    
    //First Derivative Test
    if( _parser.template get<int> ( "DerivativeTest.order" ) == 1 ){
      aol::consoleOutput( "test first Derivative" );
      VectorType testMaterial = VectorType::Random( initMaterial.size() );
      DerivativeTester<DataTypeContainer> ( totalEnergyOp, _parser.template get<RealType>( "DerivativeTest.stepSize" ), testMaterial, 1, _parser.template get<RealType> ( "DerivativeTest.tolerance" ) ).testFirstDerivative_AllDirections( );
      return;
    }
    
    cout << "5" << endl;
    SolverInfo<DataTypeContainer> solverInfo;
    switch ( _parser.template get<int> ( "MaterialOptimization.Method" ) ){
        
        // 0 - BFGS (optimize without any constraints)
        case 0 : {
                            const bool writeTimestepsBFGS = true;
                            QuasiNewtonBFGS< DataTypeContainer, ArmijoStepsizeControl > quasiNewtonSolver ( totalEnergyOp, _maxIterations, _tolerance, 10, writeTimestepsBFGS );
                            quasiNewtonSolver.solve( initMaterial, solutionMaterial );
//                             if( _parser.template get<bool> ("MaterialOptimization.printTrajectory") ){
//                                 for( int i=0; i < quasiNewtonSolver.getSolutionTrajectory().size(); i++ ){
//                                     quocHandler.plotUndeformedShellWithMaterial ( quasiNewtonSolver.getSolutionTrajectory()[i], "solTraj_" + std::to_string(i) );
//                                 }
//                             }
        }break;
        
#ifdef USE_IPOPT        
        // 10 - IPOPT (optimize without any constraints)
        case 10 : { 
                IpoptBoxConstraintFirstOrderSolver<DataTypeContainer> ipoptSolver ( totalEnergyOp, _maxIterations, _tolerance );
                ipoptSolver.solve( initMaterial, solutionMaterial );
        }break;
        
        
        // 11 - IPOPT with box constraints on material (-1 <= m <= 1)
        case 11 : { 
               IpoptBoxConstraintFirstOrderSolver<DataTypeContainer> 
               ipoptSolver ( totalEnergyOp, _maxIterations, _tolerance, _parser.template get<RealType>("MaterialOptimization.BoxConstraint_l" ), _parser.template get<RealType>("MaterialOptimization.BoxConstraint_u" ) );
              ipoptSolver.solve( initMaterial, solutionMaterial );
        }break;

        // 12 - IPOPT with volume constraint: vol = c_vol + and box constraints on material (-1 <= m <= 1)
        case 12 : {
            const RealType fixedVolHard = _parser.template get<RealType>("MaterialOptimization.VolConstraint");
            if( fixedVolHard > 1.e-16 ){
                VolumeConstraint<MatOptConfigurator> volumeConstraintOp ( OptDeformFinder.getMatOptConfigurator() );
                IpoptInterfaceNLC<DataTypeContainer> ipoptSolver ( totalEnergyOp, volumeConstraintOp, _maxIterations, _tolerance, _parser.template get<RealType>("MaterialOptimization.BoxConstraint_l"), _parser.template get<RealType>("MaterialOptimization.BoxConstraint_u" ), fixedVolHard, fixedVolHard );
                ipoptSolver.solve( initMaterial, solutionMaterial, solverInfo );
            }else{
                cout << "solution is soft material" << endl;
                for( int i=0; i<solutionMaterial.size(); i++ ) solutionMaterial[i] = -1.0;   
            }
        }break;
        
        
        // 13 - IPOPT with volume constraint: 0 \leq vol \leq c_vol and box constraints
        case 13 : {
            const RealType fixedVolHard = _parser.template get<RealType>("MaterialOptimization.VolConstraint");
            if( fixedVolHard > 1.e-16 ){
                VolumeConstraint<MatOptConfigurator> volumeConstraintOp ( OptDeformFinder.getMatOptConfigurator() );
                IpoptInterfaceNLC<DataTypeContainer> ipoptSolver ( totalEnergyOp, volumeConstraintOp, _maxIterations, _tolerance, _parser.template get<RealType>("MaterialOptimization.BoxConstraint_l"), _parser.template get<RealType>("MaterialOptimization.BoxConstraint_u" ), 0.0, fixedVolHard );
                ipoptSolver.solve( initMaterial, solutionMaterial, solverInfo );
            }else{
                cout << "solution is soft material" << endl;
                for( int i=0; i<solutionMaterial.size(); i++ ) solutionMaterial[i] = -1.0;   
            }
        }break;
#endif //USE_IPOPT
        //TODO Newton
            //cout << "start Newton" << endl;
            //RealType minTau = 1.e-7; RealType maxTau = 4.;
            //NewtonMethod<RealType,VectorType,VectorType,SparseMatrixType > NewtonSolver( gradientOp, hessianOp, initMaterial.size(), initMaterial.size(), _parser.template get<int>("MaxCGIterations"), 1.e-14, minTau, maxTau, 2 );
            //NewtonSolver.solve( initMaterial, solutionMaterial ); 
        default :
            throw std::invalid_argument( aol::strprintf ( "Wrong NonlinearOptimizationType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
            break;
    }
       
    solDisp.resize( OptDeformFinder.getSolutionDisplacement().size() ); solDisp = OptDeformFinder.getSolutionDisplacement();    
    
    //! save results
    VectorType solutionMaterialScaledTo01( solutionMaterial.size() );
    for( int i=0; i<solutionMaterial.size(); i++ ) solutionMaterialScaledTo01[i] = 0.5 * ( solutionMaterial[i] + 1.0 );
    quocHandler.plotWithMaterial( OptDeformFinder.getSolutionDisplacement(), solutionMaterialScaledTo01, "solutionDeformationForSolMaterial" );
    quocHandler.plotWithMaterial( adjointProblem.getAdjointSolution(), solutionMaterialScaledTo01, "adjointSolutionForSolMaterial" );
    quocHandler.plotUndeformedWithMaterial ( solutionMaterialScaledTo01, "solutionMaterialUndeformed" );
      
    aol::consoleOutput( "Energy of solution material : " );
    RealType solutionEnergy; VectorType solutionResidual ( solutionMaterial );
    totalEnergyOp.evaluateEnergy( solutionMaterial, solutionEnergy ); totalEnergyOp.evaluateJacobian( solutionMaterial, solutionResidual );
    totalEnergyOp.printLastEnergies(); totalEnergyOp.printLastResiduals();
    
    
    //plot parameterparser to tex-file 
//     totalEnergyOp.getLastEnergies( _lastJphysEnergy, _lastVolumeEnergy, _lastAreaEnergy );
//     _lastAreaEnergy *= matOptConf._factorInterfaceCost;
//     //TODO weighted energies 
//     //TODO volconstraint only for case 10
//     TikzPlotterMaterialOptimization<DataTypeContainer> tikzPlotter( _parser );
//     tikzPlotter.plotParameters( matOptConf, eps_area );
//     //tikzPlotter.plotSingleResultVolCompare( _lastJphysEnergy, _parser.template get<RealType>("MaterialOptimization.VolConstraint"), _lastAreaEnergy );
//     tikzPlotter.plotSingleResultVolCompare( _lastJphysEnergy, _lastVolumeEnergy, _lastAreaEnergy );
  }
  
                      
};

#endif //__QUOCMATERIALOPTIMIZATIONSOLVER_H
