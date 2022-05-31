#ifndef __QUOCOPTIMALDEFORMSOLVER_H
#define __QUOCOPTIMALDEFORMSOLVER_H

#include <general.h>
#include <loadAndSave.h>

// #include <meshWithData.h>
#include <meshWithData.h>
#include "../quocMaterialOptimizationDefines.h"
#include <linearSolver.h>
#include <quocHandler.h>

#include "quocLinElastEnergy.h"
#include "quocForces.h"

using namespace quocFE;

template <typename MatOptConfigurator>
class quocOptimalDeformSolverInterface {
protected :
  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType      ParameterParserType;
  
  const ParameterParserType &_parser;
  const MatOptConfigurator &_matOpConf;
  const ConfiguratorType & _conf;
  const QuocHandler<ConfiguratorType> & _quocHandler;

  // phasefield describing material
  mutable VectorType _pf;
  const MaskType & _mask;
  
  mutable VectorType _rhs;
  mutable VectorType _solDisp;
  
  mutable RealType _energy_Elast, _energy_Pot,  _energy_Total;

  public:

      quocOptimalDeformSolverInterface( const ParameterParserType &Parser,
                                        const MatOptConfigurator & matOpConf,
                                        const VectorType & Phasefield,
                                        const QuocHandler<ConfiguratorType> & quocHandler,
                                        const VectorType & rhs_force ) :
         _parser ( Parser ),
         _matOpConf ( matOpConf ),
         _conf ( matOpConf._conf ),
         _quocHandler ( quocHandler ),
         _pf ( Phasefield ),
         _mask ( quocHandler.getDirichletMask() ),
         _rhs ( rhs_force ),
         _solDisp ( _conf.dimDomain * _conf.getNumGlobalDofs() ){ }

    const ParameterParserType& getParser( ) const { return _parser; }
    const MatOptConfigurator & getMatOptConfigurator() const { return _matOpConf;}
    const QuocHandler<ConfiguratorType> & getQuocHandler ( ) const { return _quocHandler; }

    const VectorType & getCurrentPf ( ) const { return _pf; };
    const VectorType& getRHS ( ) const{  return _rhs;  }
    const VectorType& getSolutionDisplacement() const { return _solDisp;}

    void getLastEnergy ( RealType &energy_Elast, RealType & energy_P, RealType & energy_Total ){
     energy_Elast = _energy_Elast;
     energy_P = _energy_Pot;
     energy_Total = _energy_Total;
    }

//     void OutputLastEnergyAndResidualLagrangian ( ){
//      cout << "energy_M = "   << _energy_Membrane << endl
//           << "energy_B = "   << _energy_Bending << endl
//           << "energy_P = "   << _energy_Pot << endl;
//     }
  
};



template <typename MatOptConfigurator, typename LinElastEnergyType = LinElastHessian<MatOptConfigurator> >
class quocOptimalDeformSolverLinElast : public quocOptimalDeformSolverInterface<MatOptConfigurator> {

  typedef typename MatOptConfigurator::ConfiguratorType        ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer               DataTypeContainer;
  typedef typename ConfiguratorType::RealType                  RealType;
  typedef typename ConfiguratorType::MaskType                  MaskType;
  typedef typename ConfiguratorType::InitType                  MeshType;
  typedef typename ConfiguratorType::VectorType                VectorType;
  typedef typename ConfiguratorType::SparseMatrixType          SparseMatrixType;
  typedef typename DataTypeContainer::ParameterParserType      ParameterParserType;
  typedef typename DataTypeContainer::TripletType              TripletType;
 
protected :
  mutable SparseMatrixType _HessianLinElast;

  public:
      quocOptimalDeformSolverLinElast( const ParameterParserType &Parser,
                                       const MatOptConfigurator & matOpConf,
                                       const VectorType & Phasefield,
                           const QuocHandler<ConfiguratorType> & quocHandler,
                           const VectorType & rhs_force ) :
         quocOptimalDeformSolverInterface<MatOptConfigurator> ( Parser, matOpConf, Phasefield, quocHandler, rhs_force ),
         _HessianLinElast( matOpConf._conf.dimDomain * matOpConf._conf.getNumGlobalDofs(), matOpConf._conf.dimDomain * matOpConf._conf.getNumGlobalDofs() )
        {
          LinElastEnergyType ( this->_matOpConf, this->_pf ).assembleDirichlet( _HessianLinElast, this->_mask );
          _HessianLinElast.makeCompressed();
          solve();
        }

    const SparseMatrixType& getHessianLinElast () const { return _HessianLinElast; }
    const SparseMatrixType& getSystemMatForAdjointProblem () const { return _HessianLinElast; }
    
    void updatePhasefield( const VectorType &pf ) const {
      RealType diff =( this->_pf - pf ).squaredNorm();
      if ( diff > 1.e-15 ){
        this->_pf = pf;
        LinElastEnergyType ( this->_matOpConf, this->_pf ).assembleDirichlet( _HessianLinElast, this->_mask );
        _HessianLinElast.makeCompressed();
        solve();
      }
    }

protected:
 
  void solve() const {
    solveLinearSystem<DataTypeContainer>( _HessianLinElast, this->_solDisp, this->_rhs, static_cast<LinearSolverMethod> (this->_parser.template get<int>( "ConstraintProblem.LinearSolverMethod" ) ) );
  }
  
};

#endif
