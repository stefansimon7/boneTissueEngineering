#ifndef __QUOCMATERIALOPTIMIZATIONENERGIES_H
#define __QUOCMATERIALOPTIMIZATIONENERGIES_H

#include "../elastDeformLinElast/quocLinElastEnergy.h"
#include "quocEnergiesWithMaterial.h"
#include "quocMaterialOptimizationSolver.h"
#include "quocMaterialOptimizationComplianceEnergies.h"

using namespace quocFE;

template <typename MatOptConfigurator, typename LinElastEnergyType = LinElastHessian<MatOptConfigurator> >
class MaterialOptimizationEnergyOp : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer > 
{
  protected :
    typedef typename MatOptConfigurator::ConfiguratorType::DTContainer DataTypeContainer;
    typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
    typedef typename MatOptConfigurator::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;

    const quocOptimalDeformSolverLinElast<MatOptConfigurator,LinElastEnergyType> & _OptDeformFinder;
    const MaterialOptimizationAdjointProblem<MatOptConfigurator,LinElastEnergyType> & _adjointProblem;
    const RealType _c_compl, _c_vol, _c_interface, _epsInterfaceLength;
    
    mutable RealType _lastVolumeEnergy, _lastAreaEnergy, _lastJphysEnergy, _lastTotalEnergy;
    mutable RealType _lastVolumeResidual, _lastAreaResidual, _lastJphysResidual, _lastTotalResidual;
    
  public:
    MaterialOptimizationEnergyOp( const MatOptConfigurator &matOptConf,
                                  const quocOptimalDeformSolverLinElast<MatOptConfigurator,LinElastEnergyType> & OptDeformFinder,
                                  const MaterialOptimizationAdjointProblem<MatOptConfigurator,LinElastEnergyType> & adjointProblem ) : 
         _OptDeformFinder( OptDeformFinder ), 
         _adjointProblem ( adjointProblem ),
         _c_compl ( matOptConf._factorComplianceCost ),
         _c_vol( matOptConf._factorVolumeCost ), 
         _c_interface( matOptConf._factorInterfaceCost ),
         _epsInterfaceLength( matOptConf._epsInterfaceLength )
         { } 
   
   void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
         // compute volume term vol(v) = int 1/4(v+1)^2 = 1/4 (V^t M V + 2 V^t M 1 + 1^t M 1)
         PfOp_Volume<MatOptConfigurator> ( _OptDeformFinder.getMatOptConfigurator() ).apply( v, _lastVolumeEnergy );
         energy += _lastVolumeEnergy * _c_vol;
         // compute area term
         PfOp_ModicaMortola<MatOptConfigurator>( _OptDeformFinder.getMatOptConfigurator(), _epsInterfaceLength ).apply( v, _lastAreaEnergy); 
         energy += _lastAreaEnergy * _c_interface;
         // compute compliance
         ComplianceEnergy<MatOptConfigurator,LinElastEnergyType>( _OptDeformFinder ).apply( v, _lastJphysEnergy );
         energy += _lastJphysEnergy * _c_compl;
         _lastTotalEnergy = energy;
    }
    
    RealType getLastVolumeEnergy() const { return _lastVolumeEnergy;}
    RealType getLastAreaEnergy()   const { return _lastAreaEnergy;}
    RealType getLastJphysEnergy()  const { return _lastJphysEnergy;}
    
    void getLastEnergies( RealType &Compliance, RealType &Area, RealType &Interface ) const{ Compliance = _lastJphysEnergy; Area = _lastVolumeEnergy; Interface = _lastAreaEnergy;  }
    
    void printLastEnergies() const {
     cout << "Energy:" << endl
          << std::fixed << std::setprecision( 12 )
          << "Area       = " << getLastVolumeEnergy() << " , \t c_a Area       = " << _c_vol * getLastVolumeEnergy() << endl
          << "Interface  = " << getLastAreaEnergy()   << " , \t c_i Interface  = " << _c_interface * getLastAreaEnergy() << endl
          << "Compliance = " << getLastJphysEnergy()  << " , \t c_c Compliance = " << _c_compl * getLastJphysEnergy() << endl
          << "Total (weighted ) = " << _lastTotalEnergy << endl;
    }
   
   void evaluateJacobian( const VectorType &v, VectorType& Deriv ) const {  
         Deriv.setZero(); 
         VectorType grad ( v.size() );
         // compute volume term gradient
         PfOp_VolumeDerivative<MatOptConfigurator> ( _OptDeformFinder.getMatOptConfigurator() ).apply( v, grad );
         Deriv += _c_vol * grad;
         _lastVolumeResidual = grad.norm();
         // compute derivative of area term:
         PfOp_ModicaMortolaDerivative<MatOptConfigurator>( _OptDeformFinder.getMatOptConfigurator(), _epsInterfaceLength ).apply( v, grad );
         Deriv += _c_interface * grad;
         _lastAreaResidual = grad.norm();
         // compute gradient of Jphys w.r.t. phasefield Jphys,v 
         ComplianceEnergy_Derivative<MatOptConfigurator,LinElastEnergyType>( _OptDeformFinder, _adjointProblem ).apply( v, grad );
         Deriv += _c_compl * grad;
         _lastJphysResidual = grad.norm();
         _lastTotalResidual = Deriv.norm();
   }
   
    RealType getLastVolumeResidual() const { return _lastVolumeResidual;}
    RealType getLastAreaResidual()   const { return _lastAreaResidual;}
    RealType getLastJphysResidual()  const { return _lastJphysResidual;}
    
    void printLastResiduals() const {
     cout << "Residuals :" << endl
          << "Area       = " << getLastVolumeResidual() << endl
          << "Interface  = " << getLastAreaResidual() << endl
          << "Compliance = " << getLastJphysResidual() << endl
          << "total      = " << _lastTotalResidual << endl; 
    }
   
   void evaluateHessian( const VectorType &/*v*/, SparseMatrixType& /*Hessian*/ ) const {  
       throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
   }
   
   void evaluateTripletListHessian( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
   }
   
   void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
   }

};









#endif //__MATERIALOPTIMIZATIONENERGIES_H