#ifndef __QUOCENERGIESWITHMATERIAL_H
#define __QUOCENERGIESWITHMATERIAL_H

#include <general.h>
#include <energyDefines.h>
#include <quocIntegrator.h>
#include <quocHandler.h>
#include "../quocMaterialOptimizationDefines.h" //for approxCharFct, doubleWell

using namespace quocFE;
    
//==========================================================================================================================
// volume \int_domain chi(m)
//==========================================================================================================================
template <typename MatOptConfiguratorType>
class PfOp_Volume :
public QuocIntegrator < typename MatOptConfiguratorType::ConfiguratorType, PfOp_Volume<MatOptConfiguratorType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const MatOptConfiguratorType & _matOptConf;
  mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
  
  public:
    PfOp_Volume ( const MatOptConfiguratorType & matOptConf ) : 
     QuocIntegrator<ConfiguratorType, PfOp_Volume<MatOptConfiguratorType> > ( matOptConf._conf ),
     _matOptConf ( matOptConf ), 
     _PfPtr ( NULL ) {}
     
  ~PfOp_Volume() {
     delete _PfPtr;
  };
  
  void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocIntegrator< ConfiguratorType, PfOp_Volume<MatOptConfiguratorType>  >::assembleAdd( energy );
  }

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     return _matOptConf.approxCharFct_vol( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


// D Area(v) (hat v) = \int chi'(v) (hat v)
template <typename MatOptConfiguratorType>
class PfOp_VolumeDerivative : public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_VolumeDerivative<MatOptConfiguratorType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const MatOptConfiguratorType & _matOptConf;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_VolumeDerivative ( const MatOptConfiguratorType & matOptConf ) 
  : QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_VolumeDerivative<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _PfPtr ( NULL ) {}
     
  ~PfOp_VolumeDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_VolumeDerivative<MatOptConfiguratorType> >::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return _matOptConf.approxCharFct_vol_Derivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


// D^2 Area(v) (hat v) (dhat v) = \int chi''(v) (hat v)(dhat v)
template <typename MatOptConfiguratorType>
class PfOp_VolumeSecondDerivative : public QuocFELinWeightedMassIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_VolumeSecondDerivative<MatOptConfiguratorType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const MatOptConfiguratorType & _matOptConf;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_VolumeSecondDerivative ( const MatOptConfiguratorType & matOptConf ) 
  : QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), 
    _PfPtr ( NULL ) {}
     
  ~PfOp_VolumeSecondDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<MatOptConfiguratorType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return _matOptConf.approxCharFct_vol_SecondDerivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


// int_S chi(v)
template <typename MatOptConfigurator>
class VolumeConstraint : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  
    protected:
        const MatOptConfigurator &_matOptConf;

    public:
    VolumeConstraint ( const MatOptConfigurator &matOptConf ) : _matOptConf ( matOptConf ){ } 
      
    void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
         PfOp_Volume<MatOptConfigurator> ( _matOptConf ).apply( v, energy );
    }
    
    void evaluateJacobian( const VectorType & v, VectorType& Deriv ) const {
         PfOp_VolumeDerivative<MatOptConfigurator> ( _matOptConf ).apply( v, Deriv );
    }
    
    void evaluateHessian( const VectorType & v, SparseMatrixType& Hessian ) const {
         PfOp_VolumeSecondDerivative<MatOptConfigurator> ( _matOptConf ).apply( v, Hessian );
    }
    
    void evaluateTripletListHessian( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
      
};


template <typename MatOptConfigurator>
class VolumeConstraintPeriodicBC : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  
  protected:
        const MatOptConfigurator &_matOptConf;
        const QuocHandler<ConfiguratorType> &_quocHandler;

  public:
  VolumeConstraintPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler )  : _matOptConf ( matOptConf ), _quocHandler ( quocHandler ) { } 
      
    void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOp_Volume<MatOptConfigurator> ( _matOptConf ).apply( vCollabsedAndExtended, energy );
    }
    
    void evaluateJacobian( const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOp_VolumeDerivative<MatOptConfigurator> ( _matOptConf ).apply( vCollabsedAndExtended, Deriv );
         _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessian( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
      
};







 
//==========================================================================================================================
// barycenter_i(m) = (1/int_domain chi(m) ) \int_domain chi(m) x_i  -  c_i
//==========================================================================================================================
template <typename MatOptConfiguratorType>
class PfOp_Barycenter :
public QuocIntegrator < typename MatOptConfiguratorType::ConfiguratorType, PfOp_Barycenter<MatOptConfiguratorType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::PointType PointType;
  
  const MatOptConfiguratorType & _matOptConf;
  const int _direction; 
  const RealType _c;
  mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
  
  public:
    PfOp_Barycenter ( const MatOptConfiguratorType & matOptConf, const int direction, const RealType c ) : 
     QuocIntegrator<ConfiguratorType, PfOp_Barycenter<MatOptConfiguratorType> > ( matOptConf._conf ),
     _matOptConf ( matOptConf ), _direction ( direction ), _c ( c ),
     _PfPtr ( NULL ) {}
     
  ~PfOp_Barycenter() {
     delete _PfPtr;
  };
  
  void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocIntegrator< ConfiguratorType, PfOp_Barycenter<MatOptConfiguratorType>  >::assembleAdd( energy );
  }

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     PointType GlobalCoord; _matOptConf._conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
     return _matOptConf.approxCharFct_vol( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


template <typename MatOptConfiguratorType>
class PfOp_BarycenterDerivative : public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_BarycenterDerivative<MatOptConfiguratorType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::PointType PointType;
    
    const MatOptConfiguratorType & _matOptConf;
    const int _direction; 
    const RealType _c;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_BarycenterDerivative ( const MatOptConfiguratorType & matOptConf, const int direction, const RealType c  ) 
  : QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_BarycenterDerivative<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _direction ( direction ), _c ( c ), _PfPtr ( NULL ) {}
     
  ~PfOp_BarycenterDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_BarycenterDerivative<MatOptConfiguratorType> >::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      PointType GlobalCoord; _matOptConf._conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
      return _matOptConf.approxCharFct_vol_Derivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


// D^2 Area(v) (hat v) (dhat v) = \int chi''(v) (hat v)(dhat v)
template <typename MatOptConfiguratorType>
class PfOp_BarycenterSecondDerivative : public QuocFELinWeightedMassIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_BarycenterSecondDerivative<MatOptConfiguratorType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    typedef typename ConfiguratorType::PointType PointType;
    
    const MatOptConfiguratorType & _matOptConf;
    const int _direction; 
    const RealType _c;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_BarycenterSecondDerivative ( const MatOptConfiguratorType & matOptConf, const int direction, const RealType c  ) 
  : QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_BarycenterSecondDerivative<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _direction ( direction ), _c ( c ),  _PfPtr ( NULL ) {}
     
  ~PfOp_BarycenterSecondDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_BarycenterSecondDerivative<MatOptConfiguratorType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      PointType GlobalCoord; _matOptConf._conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
      return _matOptConf.approxCharFct_vol_SecondDerivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


// int_S chi(v)
template <typename MatOptConfigurator>
class BaryCenterConstraint : public aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
    protected:
        const MatOptConfigurator &_matOptConf;
        const PointType &_c;

    public:
    BaryCenterConstraint ( const MatOptConfigurator &matOptConf, const PointType &c ) 
    : aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  > ( c.size() ),
      _matOptConf ( matOptConf ), _c ( c ) { } 
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
         energy = 0;
         PfOp_Barycenter<MatOptConfigurator> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( v, energy );
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
         PfOp_BarycenterDerivative<MatOptConfigurator> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( v, Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & v, SparseMatrixType& Hessian ) const {
         PfOp_BarycenterSecondDerivative<MatOptConfigurator> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( v, Hessian );
    }
      
};


template <typename MatOptConfigurator>
class BaryCenterConstraintPeriodicBC : public aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
    protected:
        const MatOptConfigurator &_matOptConf;
        const QuocHandler<ConfiguratorType> & _quocHandler;
        const PointType &_c;


  public:
  BaryCenterConstraintPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler, const PointType &c ) 
  : aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  > ( c.size() ),
  _matOptConf ( matOptConf ), _quocHandler ( quocHandler ), _c ( c ) { }
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOp_Barycenter<MatOptConfigurator> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, energy );
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOp_BarycenterDerivative<MatOptConfigurator> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, Deriv );
         _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
      
};





//======================================================================================================================================
//================================= ModicaMortola and Derivative ============================================================================
//======================================================================================================================================
//! class to compute perimeter_eps(v) 
// = 1/2 \int_S eps |nabla v|^2 + 1/eps 9/16 (v^2-1)^2 d x
// Dirichlet Part
template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaPart1 : public QuocIntegrator<typename MatOptConfiguratorType::ConfiguratorType, PfOp_ModicaMortolaPart1<MatOptConfiguratorType> >{

protected:  
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    
    const MatOptConfiguratorType & _matOptConf;
    const RealType _eps;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;

public:
  PfOp_ModicaMortolaPart1 ( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) : 
     QuocIntegrator<ConfiguratorType, PfOp_ModicaMortolaPart1<MatOptConfiguratorType> >( matOptConf._conf ), 
     _matOptConf ( matOptConf ), 
     _eps ( eps_area ),
     _PfPtr ( NULL ) { }
     
  ~PfOp_ModicaMortolaPart1() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocIntegrator< ConfiguratorType, PfOp_ModicaMortolaPart1<MatOptConfiguratorType> >::assembleAdd( energy );
  }
  
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const {
      DomVecType Dv;
      _PfPtr->evaluateGradientAtQuadPoint (El, QuadPoint, Dv );
      RealType Part1 = ( Dv ).dot( Dv );                
      return 0.5 * _eps * Part1;
  }
};


template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaPart2 : public QuocIntegrator<typename MatOptConfiguratorType::ConfiguratorType, PfOp_ModicaMortolaPart2<MatOptConfiguratorType> >{

protected:  
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    
    const MatOptConfiguratorType & _matOptConf;
    
    const RealType _eps;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;

public:
  PfOp_ModicaMortolaPart2 ( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) : 
     QuocIntegrator<ConfiguratorType, PfOp_ModicaMortolaPart2<MatOptConfiguratorType> >( matOptConf._conf ), 
     _matOptConf ( matOptConf ),
     _eps ( eps_area ),
     _PfPtr ( NULL ) { }
     
  ~PfOp_ModicaMortolaPart2() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocIntegrator< ConfiguratorType, PfOp_ModicaMortolaPart2<MatOptConfiguratorType>  >::assembleAdd( energy );
  }
  
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const {
      return 0.5 * ( _matOptConf.doubleWell( _PfPtr->evaluateAtQuadPoint ( El, QuadPoint ) ) / _eps );
  }
};


template <typename MatOptConfiguratorType>
class PfOp_ModicaMortola {

  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const MatOptConfiguratorType & _matOptConf;
    const RealType _eps_area;
  
public:
  PfOp_ModicaMortola( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) : 
       _matOptConf ( matOptConf ), _eps_area( eps_area ) { }
  
  void apply( const VectorType& v, RealType& energy ) const {
    RealType temp = 0;
    PfOp_ModicaMortolaPart1<MatOptConfiguratorType> ( _matOptConf, _eps_area ).apply( v, temp );
    energy = temp;
    PfOp_ModicaMortolaPart2<MatOptConfiguratorType> ( _matOptConf, _eps_area ).apply( v, temp );
    energy += temp;
  }
  
};



//! class to compute \int \sqrt(eps nabla m nabla b_i
template<typename MatOptConfiguratorType>
class PfOp_ModicaMortolaDerivativePart1 :
public QuocFENonlinDiffOpIntegrator <typename MatOptConfiguratorType::ConfiguratorType, PfOp_ModicaMortolaDerivativePart1 <MatOptConfiguratorType> >
{      
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
  
    const MatOptConfiguratorType & _matOptConf;
    const RealType _eps;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;

  public:
    PfOp_ModicaMortolaDerivativePart1 ( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) :
     QuocFENonlinDiffOpIntegrator <ConfiguratorType, PfOp_ModicaMortolaDerivativePart1 <MatOptConfiguratorType> > ( matOptConf._conf ),
      _matOptConf ( matOptConf ), _eps ( eps_area ),   _PfPtr ( NULL ) { }
      
     ~PfOp_ModicaMortolaDerivativePart1() {
           delete _PfPtr;
     };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFENonlinDiffOpIntegrator <ConfiguratorType, PfOp_ModicaMortolaDerivativePart1 <MatOptConfiguratorType> > ::assembleAdd( Deriv );
  }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, DomVecType &NL) const {     
      _PfPtr->evaluateGradientAtQuadPoint( El, QuadPoint, NL );
      NL *= _eps; 
      
    }
}; 



//! doubleWell part
template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaDerivativePart2 : 
public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_ModicaMortolaDerivativePart2<MatOptConfiguratorType> > {

  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    const MatOptConfiguratorType & _matOptConf;
    
    const RealType _eps;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;

  public:
  PfOp_ModicaMortolaDerivativePart2 ( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) 
  : QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_ModicaMortolaDerivativePart2<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _eps ( eps_area ), _PfPtr ( NULL ) { }
      
     ~PfOp_ModicaMortolaDerivativePart2() {
           delete _PfPtr;
     };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFENonlinOpIntegrator <ConfiguratorType, PfOp_ModicaMortolaDerivativePart2 <MatOptConfiguratorType> > ::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El,  int QuadPoint ) const {   
      return 0.5 * _matOptConf.doubleWellDerivative( _PfPtr->evaluateAtQuadPoint ( El, QuadPoint ) ) / _eps;
  }
};


//! ModicaMortola Derivative
template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaDerivative {

  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const MatOptConfiguratorType & _matOptConf;
    const RealType _eps_area;
  
public:
  PfOp_ModicaMortolaDerivative( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) : 
       _matOptConf ( matOptConf ),  _eps_area( eps_area ) { }
  
  void apply( const VectorType& v, VectorType& Deriv ) const {
    VectorType temp ( Deriv.size() );
    PfOp_ModicaMortolaDerivativePart1<MatOptConfiguratorType> ( _matOptConf, _eps_area ).apply( v, temp );
    Deriv = temp;
    PfOp_ModicaMortolaDerivativePart2<MatOptConfiguratorType> ( _matOptConf, _eps_area ).apply( v, temp );
    Deriv += temp;
  }
  
};





/*


// Dirichlet-Part
template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaSecondDerivativePart1 
: public UnitTriangleFELinAsymMatrixWeightedStiffIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart1<MatOptConfiguratorType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::Matrix22 Matrix22;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const MatOptConfiguratorType & _matOptConf;
    
    const RealType _eps_area;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_ModicaMortolaSecondDerivativePart1 ( const MatOptConfiguratorType & matOptConf, const RealType eps_area ) 
  : UnitTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart1<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), 
    _eps_area ( eps_area ),
    _PfPtr ( NULL ) {}
     
  ~PfOp_ModicaMortolaSecondDerivativePart1() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      UnitTriangleFELinAsymMatrixWeightedStiffIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart1<MatOptConfiguratorType> >::assemble( Hessian );
  }
  
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, Matrix22 &Matrix ) const {
    Matrix *= * _eps_area;
  }
  
};

// Double-Well-Part
template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaSecondDerivativePart2 : 
public QuocFELinWeightedMassIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart2<MatOptConfiguratorType> > {
// ,public aol::Op< typename MatOptConfiguratorType::ConfiguratorType::VectorType, typename MatOptConfiguratorType::ConfiguratorType::SparseMatrixType > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const MatOptConfiguratorType & _matOptConf;
    
    const RealType _eps_area;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_ModicaMortolaSecondDerivativePart2 ( const MatOptConfiguratorType & matOptConf,  const RealType eps_area ) 
  : QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart2<MatOptConfiguratorType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), 
    
    _eps_area ( eps_area ),
    _PfPtr ( NULL ) {}
     
  ~PfOp_ModicaMortolaSecondDerivativePart2() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_ModicaMortolaSecondDerivativePart2<MatOptConfiguratorType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return 0.5 * _matOptConf.doubleWellSecondDerivative ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) / _eps_area;
  }
};

//! ModicaMortola Second Derivative
template <typename MatOptConfiguratorType>
class PfOp_ModicaMortolaSecondDerivative {

  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const MatOptConfiguratorType & _matOptConf;
    
    const RealType _eps_area;
  
public:
  PfOp_ModicaMortolaSecondDerivative( const MatOptConfiguratorType & matOptConf,
                                             
                                             const RealType eps_area ) : 
       _matOptConf ( matOptConf ),         
       _eps_area( eps_area ) {}
  
  void apply( const VectorType& v, SparseMatrixType& Hessian ) const {
    SparseMatrixType temp ( Hessian.rows(), Hessian.cols() );
    PfOp_ModicaMortolaSecondDerivativePart1<MatOptConfiguratorType> ( _matOptConf,_eps_area ).apply( v, temp );
    Hessian = temp;
    PfOp_ModicaMortolaSecondDerivativePart2<MatOptConfiguratorType> ( _matOptConf, _eps_area ).apply( v, temp );
    Hessian += temp;
  }
  
};*/


#endif