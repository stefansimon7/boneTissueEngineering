#ifndef __QUOCFORCES_H
#define __QUOCFORCES_H

#include <quocIntegrator.h>

using namespace quocFE;

template< typename ConfiguratorType >
class quocCreateConstantLoad :
public QuocFENonlinVectorOpIntegrator <ConfiguratorType, quocCreateConstantLoad<ConfiguratorType> >
{      
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::PointType PointType;
  
    const PointType &_f;

  public:
    quocCreateConstantLoad ( const ConfiguratorType &conf,
                         const PointType &f ) :
     QuocFENonlinVectorOpIntegrator <ConfiguratorType, quocCreateConstantLoad <ConfiguratorType> > (conf),
     _f ( f ){ }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint, PointType &NL) const {   
      NL = _f;
    }
};

template< typename ConfiguratorType >
class quocCreateBoundaryLoad :
public QuocFENonlinVectorOpIntegratorWithRefCoord <ConfiguratorType, quocCreateBoundaryLoad<ConfiguratorType> >
{      
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::PointType PointType;
  
    const ConfiguratorType &_conf;
    const PointType &_f;
    const int _direction; const RealType _bound;

  public:
    quocCreateBoundaryLoad ( const ConfiguratorType &conf,
                         const PointType &f,
                         const int direction, const RealType bound
                       ) :
     QuocFENonlinVectorOpIntegratorWithRefCoord <ConfiguratorType, quocCreateBoundaryLoad <ConfiguratorType> > (conf),
     _conf(conf),
     _f ( f ),
     _direction(direction), _bound(bound) { }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int /*QuadPoint*/, const PointType &RefCoord, PointType &NL) const {  
        PointType GlobalCoord;
        _conf.getGlobalCoords( El, RefCoord, GlobalCoord );
        if( GlobalCoord[_direction] > _bound ) NL = _f;
        else NL.setZero();
    }
};


template< typename ConfiguratorType >
class quocCreatePartialBoundaryLoad :
public QuocFENonlinVectorOpIntegratorWithRefCoord <ConfiguratorType, quocCreatePartialBoundaryLoad<ConfiguratorType> >{      
  
  protected:
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::PointType PointType;
  
    const ConfiguratorType &_conf;
    const PointType &_f;

  public:
    quocCreatePartialBoundaryLoad ( const ConfiguratorType &conf, const PointType &f ) :
     QuocFENonlinVectorOpIntegratorWithRefCoord <ConfiguratorType, quocCreatePartialBoundaryLoad <ConfiguratorType> > (conf),
     _conf(conf),
     _f ( f ){ }
      
    void getNonlinearity ( const typename ConfiguratorType::ElementType &El, int /*QuadPoint*/, const PointType &RefCoord, PointType &NL) const {  
        PointType GlobalCoord;
        _conf.getGlobalCoords( El, RefCoord, GlobalCoord );
        if( (GlobalCoord[0] > 1.95) && (GlobalCoord[1] < 0.55) && (GlobalCoord[1] > 0.45) ) NL = _f;
        else NL.setZero();
    }
};


template <typename ConfiguratorType>
void quocCreateForce ( const typename ConfiguratorType::DTContainer::ParameterParserType & parser,
                       const ConfiguratorType & conf,
                       typename ConfiguratorType::VectorType & rhs_force,
                       const typename ConfiguratorType::MaskType &boundaryMask,
                       const bool collapseBoundaryValues = true ) {
  
    typedef typename ConfiguratorType::PointType PointType;
    typedef typename ConfiguratorType::RealType RealType;
    
    int numGlobalDofs = conf.getNumGlobalDofs();
    rhs_force.setZero();
    switch( parser.template get<int>( "Force.Type" ) ){
      case 0:{  
        PointType f; parser.template getFixSizeVector<RealType,PointType> ( "Force.constantLoad", f );
        quocCreateConstantLoad<ConfiguratorType>( conf, f ).assembleAdd( rhs_force );
      } break;
      case 10:{  
        PointType f; parser.template getFixSizeVector<RealType,PointType> ( "Force.constantLoad", f );
        quocCreateBoundaryLoad<ConfiguratorType>( conf, f, parser.template get<int>("Force.direction"), parser.template get<double>("Force.bound") ).assembleAdd( rhs_force );
      } break;
      case 11:{  
        PointType f; parser.template getFixSizeVector<RealType,PointType> ( "Force.constantLoad", f );
        quocCreatePartialBoundaryLoad<ConfiguratorType>( conf, f ).assembleAdd( rhs_force );
      } break;
      default:
          throw std::invalid_argument( aol::strprintf ( "Unknown ForceType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
    }

    //bc rhs_force 
    if( collapseBoundaryValues ){
        for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<ConfiguratorType::dimDomain; ++comp )
              rhs_force[i + comp * numGlobalDofs] = 0.0;
          }
        } 
    }
}
#endif
