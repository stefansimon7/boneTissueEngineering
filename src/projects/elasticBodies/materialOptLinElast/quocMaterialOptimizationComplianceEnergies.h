#ifndef __QUOCMATERIALOPTIMIZATIONCOMPLIANCEENERGIES_H
#define __QUOCMATERIALOPTIMIZATIONCOMPLIANCEENERGIES_H

#include <quocHandler.h>
#include "../quocMaterialOptimizationDefines.h"
#include "../elastDeformLinElast/quocOptimalDeformSolverLinElast.h"

using namespace quocFE;

//======================================================================================================================================
//================================= Interfaces  ===============================================================================
//======================================================================================================================================


// template<typename MatOptConfiguratorType, typename LinElastEnergyType>
// class ComplianceEnergy { };



//======================================================================================================================================
//================================= potential energy Jphys(m) = \int_M f \cdot u[m] dx  ===============================================
//======================================================================================================================================

template<typename MatOptConfiguratorType, typename LinElastEnergyType = LinElastHessian<MatOptConfiguratorType> >
class ComplianceEnergy {
  
  typedef typename MatOptConfiguratorType::RealType RealType;
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::VectorType VectorType;

protected:
  const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & _OptDeformFinder;
  
public:
  ComplianceEnergy( const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & OptDeformFinder  ) : _OptDeformFinder ( OptDeformFinder ){}

  void apply( const VectorType & v, RealType& energy ) const { 
     _OptDeformFinder.updatePhasefield( v );
      energy = ( _OptDeformFinder.getRHS() ).dot( _OptDeformFinder.getSolutionDisplacement() );
  }
  
};


//======================================================================================================================================
//================================= Derivative of Jphys  ===============================================================================
//======================================================================================================================================


//======================================================================================================================================
//================================= D_m \hat Jphys(m,u)  ===============================================================================
//======================================================================================================================================
// template <typename MatOptConfiguratorType, typename LinElastEnergyType >
// class ComplianceEnergy_DerivativeInM {};


template <typename MatOptConfiguratorType, typename LinElastEnergyType = LinElastHessian<MatOptConfiguratorType> >
class ComplianceEnergy_DerivativeInM {
   
protected: 
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::VectorType VectorType;

public:
    ComplianceEnergy_DerivativeInM ( const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & /*OptDeformFinder*/ ) {}
    void apply ( const VectorType & /*Pf*/, VectorType &Dest ) const { Dest.setZero(); }
};

//======================================================================================================================================
//================================= D_u \hat Jphys(m,u) (used for RHS in adjoint Problem) ==============================================
//======================================================================================================================================


template < typename MatOptConfiguratorType, typename LinElastEnergyType = LinElastHessian<MatOptConfiguratorType> >
class ComplianceEnergy_DerivativeInDisp {

  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & _OptDeformFinder;
  
public:
  ComplianceEnergy_DerivativeInDisp( const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & OptDeformFinder  ) :
  _OptDeformFinder ( OptDeformFinder ){}
  
  template <typename VectorTypeRef = VectorType &>
  void assemble( VectorTypeRef Deriv ) const { Deriv = _OptDeformFinder.getRHS(); }
};


//======================================================================================================================================
//================================= Adjoint Problem ==============================================
//======================================================================================================================================

template <typename MatOptConfiguratorType, typename LinElastEnergyType = LinElastHessian<MatOptConfiguratorType> >
class MaterialOptimizationAdjointProblem {
protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> &_OptDeformFinder;
    mutable VectorType _adjointSolution;

public:
    MaterialOptimizationAdjointProblem ( const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> &OptDeformFinder ) : 
     _OptDeformFinder ( OptDeformFinder ),
     _adjointSolution ( -1. * _OptDeformFinder.getSolutionDisplacement() ) { }

    void update ( ) const { _adjointSolution = -1. * _OptDeformFinder.getSolutionDisplacement(); }
    const VectorType& getAdjointSolution () const { return _adjointSolution; }
};

//======================================================================================================================================
//================================= D^2 E ( m, u(m) ) AdjointSolution  ===============================================================================
//======================================================================================================================================

//! D^2_(m,u) E(m, u(m)) (hat m)(P)
template<typename MatOptConfiguratorType, typename LinElastEnergyType = LinElastHessian<MatOptConfiguratorType> >
class ComplianceEnergy_DerivativePartAdjoint
: public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, ComplianceEnergy_DerivativePartAdjoint<MatOptConfiguratorType, LinElastEnergyType> > {

protected: 
  
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
//   typedef typename LinElastEnergyType::EvaluationHelper EvaluationHelper; //TODO 
  
  const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> &_OptDeformFinder;
  const MaterialOptimizationAdjointProblem<MatOptConfiguratorType,LinElastEnergyType> &_adjointProblem;
  const Material<RealType> & _HardMaterial;
  const RealType _factorVoidMaterial;
//   const EvaluationHelper _LinElastHelper; //TODO
  mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
  mutable QuocDiscreteVectorFunctionDefault<ConfiguratorType> *_dispPtr, *_adjointPtr;
  
public:
  ComplianceEnergy_DerivativePartAdjoint ( const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & OptDeformFinder, 
                                           const MaterialOptimizationAdjointProblem<MatOptConfiguratorType,LinElastEnergyType> &adjointProblem ) 
 : QuocFENonlinOpIntegrator<ConfiguratorType, ComplianceEnergy_DerivativePartAdjoint<MatOptConfiguratorType,LinElastEnergyType> > ( OptDeformFinder.getMatOptConfigurator()._conf ), 
  _OptDeformFinder ( OptDeformFinder ), _adjointProblem ( adjointProblem ),
   _HardMaterial ( OptDeformFinder.getMatOptConfigurator()._materialInfo._HardMaterial ),
   _factorVoidMaterial( OptDeformFinder.getMatOptConfigurator()._factorVoidMaterial ),
//   _LinElastHelper ( OptDeformFinder.getMatOptConfigurator(), OptDeformFinder.getShellHandler().getChartToUndeformedShell_Cache() ),
  _PfPtr ( NULL ), _dispPtr( NULL ), _adjointPtr( NULL ) {}
      
   ~ComplianceEnergy_DerivativePartAdjoint() {
      delete _PfPtr;
      delete _dispPtr;
      delete _adjointPtr;
   };

  void apply ( const VectorType &Pf, VectorType &Dest ) const {
    Dest.setZero();
    _OptDeformFinder.updatePhasefield( Pf ); // This computes opt disp u(Pf) and by reference this updated in adjoint problem
    _adjointProblem.update( ); // this computes solution of adjoint problem
    delete _dispPtr;
    _dispPtr = new QuocDiscreteVectorFunctionDefault<ConfiguratorType> ( _OptDeformFinder.getMatOptConfigurator()._conf, _OptDeformFinder.getSolutionDisplacement() );
    delete _adjointPtr;
    _adjointPtr = new QuocDiscreteVectorFunctionDefault<ConfiguratorType> ( _OptDeformFinder.getMatOptConfigurator()._conf, _adjointProblem.getAdjointSolution() );
    delete _PfPtr;
    _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _OptDeformFinder.getMatOptConfigurator()._conf, Pf );
    QuocFENonlinOpIntegrator< ConfiguratorType, ComplianceEnergy_DerivativePartAdjoint<MatOptConfiguratorType,LinElastEnergyType> >::assembleAdd( Dest );
  }

  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      //TODO so far only for LinElastHessian
//     RealType NL = _LinElastHelper.applyLocalMatrix_Gradient( El, QuadPoint, _OptDeformFinder.getSolutionDisplacement(), _adjointProblem.getAdjointSolution(), _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
    RealType div_disp, div_adj;
    div_disp = _dispPtr->evaluateDivergenceAtQuadPoint( El, QuadPoint );
    div_adj = _adjointPtr->evaluateDivergenceAtQuadPoint( El, QuadPoint );
    DerivativeVectorValuedType GradSym_disp, GradSym_adj;
    _dispPtr->evaluateSymmetrizedGradientAtQuadPoint( El, QuadPoint, GradSym_disp );
    _adjointPtr->evaluateSymmetrizedGradientAtQuadPoint( El, QuadPoint, GradSym_adj );
    RealType Dchi = _OptDeformFinder.getMatOptConfigurator().approxCharFct_material_Derivative( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
    RealType materialfactor = Dchi * (1. - _factorVoidMaterial);
    return materialfactor * ( _HardMaterial.getLambda() * div_disp * div_adj + 2. * _HardMaterial.getMu() * aol::ddProd<RealType,DerivativeVectorValuedType>(GradSym_disp, GradSym_adj) );
      
  }

};


//======================================================================================================================================
//================================= D_m Jphys(m)  ===============================================================================
//======================================================================================================================================


template<typename MatOptConfiguratorType, typename LinElastEnergyType = LinElastHessian<MatOptConfiguratorType> >
class ComplianceEnergy_Derivative {
    
protected: 
  
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
    
  const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> &_OptDeformFinder;
  const MaterialOptimizationAdjointProblem<MatOptConfiguratorType,LinElastEnergyType> &_adjointProblem;
  
public:
  
  ComplianceEnergy_Derivative ( const quocOptimalDeformSolverLinElast<MatOptConfiguratorType,LinElastEnergyType> & OptDeformFinder, 
                                const MaterialOptimizationAdjointProblem<MatOptConfiguratorType,LinElastEnergyType> &adjointProblem ) :
        _OptDeformFinder ( OptDeformFinder ), _adjointProblem ( adjointProblem ) {}
  
  void apply( const VectorType &Pf, VectorType &Dest ) const {
      
      VectorType tmp1 ( Dest.size() ), tmp2 ( Dest.size() );
      
      ComplianceEnergy_DerivativeInM<MatOptConfiguratorType, LinElastEnergyType> ( _OptDeformFinder ).apply( Pf, tmp1 );
      ComplianceEnergy_DerivativePartAdjoint<MatOptConfiguratorType, LinElastEnergyType> ( _OptDeformFinder, _adjointProblem ).apply( Pf, tmp2 );
      
      Dest = tmp1 + tmp2;
  }
    
    
};


#endif //__MATERIALOPTIMIZATIONCOMPLIANCEENERGIES_H