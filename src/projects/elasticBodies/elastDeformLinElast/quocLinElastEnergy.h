#ifndef __QUOCLINELASTENERGY_H
#define __QUOCLINELASTENERGY_H

#include <quocIntegrator.h>
#include "../quocMaterialOptimizationDefines.h"

//TODO int domain replaced by dim from mesh 

using namespace quocFE;

template<typename MatOpType>
class LinElastEnergy :
public QuocIntegrator<typename MatOpType::ConfiguratorType, LinElastEnergy<MatOpType> >
{
  protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType        DerivativeVectorValuedType;
    
    const MatOpType &_matOpConf;
    const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
    const QuocDiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::dimDomain> _displacement;
    const Material<RealType> & _HardMaterial;
  const RealType _factorVoidMaterial;

    
  public:
    LinElastEnergy ( const MatOpType &matOpConf, const VectorType &material, const VectorType &disp ) : 
     QuocIntegrator <ConfiguratorType, LinElastEnergy<MatOpType> > ( matOpConf._conf ), 
     _matOpConf(matOpConf), 
     _pf( matOpConf._conf, material ),
     _displacement ( matOpConf._conf, disp ),
     _HardMaterial ( matOpConf._materialInfo._HardMaterial ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) { }
       
    RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int quadPoint ) const{
       RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
       RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
       RealType divergence = _displacement.evaluateDivergenceAtQuadPoint(El, quadPoint );
       DerivativeVectorValuedType symGrad; _displacement.evaluateSymmetrizedGradientAtQuadPoint(El,quadPoint,symGrad);
       RealType aux = 0.0;
       aux +=  _HardMaterial.getLambda() * aol::Sqr( divergence );
       aux += 2. * _HardMaterial.getMu() * symGrad.squaredNorm();
       return materialfactor * aux;
    }
};


//! The corresponding matrix assembly yields 
//! \f$ \left(\int_\Omega  lambda div(\phi_i) div(\phi_j) + 2 mu eps(u):eps(v) dx\right)_{ij} \f$ 
//! for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename MatOpType>
class LinElastHessian :
      public QuocBlockMatrixValuedIntegratorBase< typename MatOpType::ConfiguratorType, LinElastHessian<MatOpType> > {
protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::VectorType VectorType;
  const MatOpType &_matOpConf;
  const ConfiguratorType &_config;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
  const Material<RealType> & _HardMaterial;
  const RealType _factorVoidMaterial;
//   const Material<RealType> & _SoftMaterial;

public:
  LinElastHessian ( const MatOpType &matOpConf, const VectorType &material ) : 
   QuocBlockMatrixValuedIntegratorBase< ConfiguratorType, LinElastHessian<MatOpType> > ( matOpConf._conf ),
   _matOpConf( matOpConf ),
   _config ( matOpConf._conf ), 
   _pf( _config, material ),
   _HardMaterial ( matOpConf._materialInfo._HardMaterial ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) {}

  void prepareLocalMatrix( const typename ConfiguratorType::ElementType &El, LocalMatrixType (&localMatrix)[ConfiguratorType::dimDomain][ConfiguratorType::dimDomain]) const {
      
      for (int argComp = 0; argComp < ConfiguratorType::dimDomain; ++argComp)
          for (int destComp = 0; destComp < ConfiguratorType::dimDomain; ++destComp)
              localMatrix[argComp][destComp].setZero();
            
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
      const int numDofs = _config.getNumLocalDofs ( El );
      RealType nonlinearity = 0.;    
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {

          RealType chi = _matOpConf.approxCharFct_material ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int i = 0; i < numDofs; ++i ) {
              const DomVecType& grad_b_i = bfs.evaluateGradient( i, quadPoint );
              for ( int j = 0; j < numDofs; ++j ) {
                  const DomVecType& grad_b_j = bfs.evaluateGradient( j, quadPoint );
                  for(int argComp=0; argComp < ConfiguratorType::dimDomain; ++argComp){
                      localMatrix[argComp][argComp](j,i) +=  materialfactor * _HardMaterial.getMu() * grad_b_i.dot(grad_b_j) * bfs.getWeight ( quadPoint );
                      for(int destComp=0; destComp < ConfiguratorType::dimDomain; ++destComp ){
                        localMatrix[argComp][destComp](j,i) += materialfactor * _HardMaterial.getLambda() * grad_b_i[argComp] * grad_b_j[destComp] * bfs.getWeight ( quadPoint );
                        localMatrix[argComp][destComp](j,i) += materialfactor * _HardMaterial.getMu() * grad_b_i[destComp] * grad_b_j[argComp] * bfs.getWeight ( quadPoint );
                      }
                  }
              }
          }
      }
  }

};


#endif