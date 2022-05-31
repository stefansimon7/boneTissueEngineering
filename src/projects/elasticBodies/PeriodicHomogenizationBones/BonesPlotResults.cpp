#include <quocConfigurators.h>
#include <quocHandler.h>
#include <quocDiscreteFunction.h>
#include <quocDefines.h>
#include <meshWithData.h>
#include <VTKInterfaceToSurface.h>
#include "BonesPlotter.h"
#include "BonesInterface.h"
#include "BonesEnergiesWithMaterial.h"

using namespace quocFE;
using namespace shapeOptBonePolymerPeriodicHomogenization;

// #define _USEQUOC2D
#define _USEQUOC3D

#ifdef _USEQUOC2D
typedef Quoc2DDataTypeContainer                                                                                         DataTypeContainer;
typedef QuocConfigurator2D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

#ifdef _USEQUOC3D
typedef Quoc3DDataTypeContainer                                                                                         DataTypeContainer;
typedef QuocConfigurator3D<DataTypeContainer>                                                                           ConfiguratorType;
#endif

typedef typename DataTypeContainer::RealType                                                                            RealType;
typedef typename DataTypeContainer::IntVecType                                                                          IntVecType;
typedef typename DataTypeContainer::VectorType                                                                          VectorType;
typedef typename DataTypeContainer::ParameterParserType                                                                 ParameterParserType;
typedef typename DataTypeContainer::PointType                                                                           PointType;
typedef typename ConfiguratorType::SparseMatrixType                                                                     SparseMatrixType;
typedef typename ConfiguratorType::MaskType                                                                             MaskType;
typedef typename ConfiguratorType::InitType                                                                             MeshType;
typedef QuocMaterialOptimizationConfiguratorBones< ConfiguratorType >                                                   MatOptConfigurator;


//--------------------------------------------------------
int main(int argc, char ** argv) {
   
    cout << "Usage of programm BonesPlotResults:" << endl
         << "  1 - [saveDirectory]" << endl
         << "  2 - [initialization (1) or solution (2) or both (3)] (optional, default initial (1) )" << endl 
         << "  3 - [plotBlock - int numBlocksPerDirection] (optional, default false (<=1) )" << endl
         << "  4 - [plotInterface] (optional, default false)" << endl
         << endl;
  
    const string saveDirectory = argv[1];
 
    int InitialOrSolution = 1;
    if( argc > 2 ){
        InitialOrSolution = stoi( argv[2] );
    }
    
    bool plotBlock = false; int numBlocksPerDirection = 1;
    if( argc > 3 ){
        numBlocksPerDirection = stoi( argv[3] );
        if( numBlocksPerDirection > 1 ) plotBlock = true;
    }
    
    bool plotInterface = false;
    if( argc > 4 ){
        int plotInterfaceInt = stoi( argv[4] );
        if( plotInterfaceInt == 1 ) plotInterface = true;
    }
    
    bool computeHomogenizedTensor = true;
    
        
        ParameterParserType parser( aol::strprintf( "%s/ParameterParser.ini", argv[1] ).c_str() );
        parser.set ( "saving.saveDirectory", argv[1] );
    
        //=====================================================================================================================
        //Plot vtk files to png images
        //=====================================================================================================================
        aol::consoleOutput(  "plot vtk files to images" );
        const int numLoads = parser.template get<int> ("AffineDisp.numLoads");
        VTKPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> vtkPlotter( parser, numLoads );
        switch( InitialOrSolution ){
            case 1:{
                vtkPlotter.plotAll("InitMaterial");  
            }break;
            case 2:{
                vtkPlotter.plotAll("SolMaterial"); 
            }break;
            case 3:{
                vtkPlotter.plotAll();
            }break;
            default :
                throw std::invalid_argument( aol::strprintf ( "Wrong type. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
        }
        
        
        //=====================================================================================================================
        //initialize config, quocHandler
        //=====================================================================================================================
        IntVecType numDofVec; parser.template getFixSizeVector<int,IntVecType> ("InputMesh.NumDofVec", numDofVec );
        PointType lengthVec; parser.template getFixSizeVector<RealType, PointType> ("InputMesh.LengthVec", lengthVec );
        MeshType mesh ( numDofVec, lengthVec );
        ConfiguratorType conf ( mesh );
        QuocHandler<ConfiguratorType> quocHandler( parser, conf );
        
        VectorType material ( mesh.getNumVertices() );
        aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s/SolMaterial.txt", saveDirectory.c_str()  ) ); 
        quocHandler.extendVectorPeriodically( material );
        QuocDiscreteFunctionDefault<ConfiguratorType> discreteFctSingleCell ( conf, material ); //note: here solMaterialBlock is actually OldLevel
        
        const int numDofsOffsetOutside = 3;
        PointType offsetOutside; for( int i=0; i<offsetOutside.size(); ++i ) offsetOutside[i] = numDofsOffsetOutside * mesh.getMeshSize( i );
        
        //! ========================================================================================
        //! plot material (B or P) around single micro cell
        //! ========================================================================================
            aol::consoleOutput(  "plot single cell with material outside" );
            IntVecType numDofVecOutside;
            PointType lengthVecOutside;
            for( int i=0; i<numDofVec.size(); ++i ){
                numDofVecOutside[i] = numDofVec[i] + 2 * numDofsOffsetOutside;
                lengthVecOutside[i] = lengthVec[i] + 2 * numDofsOffsetOutside * mesh.getMeshSize( i ); 
            }
            MeshType meshOutside ( numDofVecOutside, lengthVecOutside );
            VectorType materialOutsideBone ( meshOutside.getNumVertices() );
            for( int i=0; i<materialOutsideBone.size(); ++i ) materialOutsideBone[i] = 1.0;
            VectorType materialOutsidePolymer ( meshOutside.getNumVertices() );
            for( int i=0; i<materialOutsidePolymer.size(); ++i ) materialOutsidePolymer[i] = -1.0;
            for( int nodeIdxOffset=0; nodeIdxOffset < meshOutside.getNumVertices(); nodeIdxOffset++ ){
                const PointType& GlobalCoordsOffset = meshOutside.getVertex ( nodeIdxOffset );
               
                bool inside = true;
                for( int i=0; i<GlobalCoordsOffset.size(); ++i ){
                  if( GlobalCoordsOffset[i] < offsetOutside[i] ) inside = false;
                  if( GlobalCoordsOffset[i] > lengthVec[i] + offsetOutside[i] ) inside = false;
                }
                
                if( inside ) {
                    PointType GlobalCoordsSingleCell;
                    for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                        GlobalCoordsSingleCell[i] = GlobalCoordsOffset[i] - numDofsOffsetOutside * mesh.getMeshSize(i);
                    }
                    int elementNumberSingleCell; PointType LocalCoordSingleCell;
                    conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                    materialOutsideBone[nodeIdxOffset] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                    materialOutsidePolymer[nodeIdxOffset] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                }
            }
            
            MeshWithData<MeshType> meshSaverOutsideBone ( meshOutside );
            meshSaverOutsideBone.addScalarData ( materialOutsideBone, "SolMaterialOutsideBone", VERTEX_DATA );
            meshSaverOutsideBone.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_OutsideBone.vtk", saveDirectory.c_str() ), offsetOutside );
            
            MeshWithData<MeshType> meshSaverOutsidePolymer ( meshOutside );
            meshSaverOutsidePolymer.addScalarData ( materialOutsidePolymer, "SolMaterialOutsidePolymer", VERTEX_DATA );
            meshSaverOutsidePolymer.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_OutsidePolymer.vtk", saveDirectory.c_str() ), offsetOutside );
        
        
        
        //=====================================================================================================================
        //Plot block of mesh with material 
        //=====================================================================================================================
        if( plotBlock ){
            aol::consoleOutput( "plot block of cells" );
            //! ========================================================================================
            //! plot block
            //! ========================================================================================
            IntVecType numDofVecBlock;
            PointType lengthVecBlock;
            for( int i=0; i<numDofVec.size(); ++i ){
                int oldSize = numDofVec[i];
                numDofVecBlock[i] = static_cast<int> ( numBlocksPerDirection * (oldSize - 1) + 1 );
                lengthVecBlock[i] = static_cast<RealType> ( numBlocksPerDirection ) * lengthVec[i]; 
            }
            MeshType meshBlock ( numDofVecBlock, lengthVecBlock );
            VectorType materialBlock ( meshBlock.getNumVertices() );
            for( int nodeIdxBlock=0; nodeIdxBlock < meshBlock.getNumVertices(); nodeIdxBlock++ ){
                const PointType& GlobalCoordsBlock = meshBlock.getVertex ( nodeIdxBlock );
                PointType GlobalCoordsSingleCell;
                for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                    int multDirection = static_cast<int> ( GlobalCoordsBlock[i] / lengthVec[i] );
                    GlobalCoordsSingleCell[i] = GlobalCoordsBlock[i] - multDirection * lengthVec[i];
                }
                int elementNumberSingleCell; PointType LocalCoordSingleCell;
                conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                materialBlock[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
            }
            MeshWithData<MeshType> meshSaver ( meshBlock );
            meshSaver.addScalarData ( materialBlock, "SolMaterialBlock", VERTEX_DATA );
            meshSaver.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_Block.vtk", saveDirectory.c_str() ) );
            
            //! ========================================================================================
            //! plot material (B or P) around block of micro cells
            //! ========================================================================================
            IntVecType numDofVecBlockOutside;
            PointType lengthVecBlockOutside;
            for( int i=0; i<numDofVec.size(); ++i ){
                numDofVecBlockOutside[i] = numDofVecBlock[i] + 2 * numDofsOffsetOutside;
                lengthVecBlockOutside[i] = lengthVecBlock[i] + 2 * numDofsOffsetOutside * mesh.getMeshSize( i ); 
            }
            MeshType meshBlockOutside ( numDofVecBlockOutside, lengthVecBlockOutside );
            VectorType materialBlockOutsideBone ( meshBlockOutside.getNumVertices() );
            for( int i=0; i<materialBlockOutsideBone.size(); ++i ) materialBlockOutsideBone[i] = 1.0;
            VectorType materialBlockOutsidePolymer ( meshBlockOutside.getNumVertices() );
            for( int i=0; i<materialBlockOutsidePolymer.size(); ++i ) materialBlockOutsidePolymer[i] = -1.0;
            for( int nodeIdxBlock=0; nodeIdxBlock < meshBlockOutside.getNumVertices(); nodeIdxBlock++ ){
                const PointType& GlobalCoordsBlockOffset = meshBlockOutside.getVertex ( nodeIdxBlock );
                
                bool inside = true;
                for( int i=0; i<GlobalCoordsBlockOffset.size(); ++i ){
                  if( GlobalCoordsBlockOffset[i] < offsetOutside[i] ) inside = false;
                  if( GlobalCoordsBlockOffset[i] > lengthVecBlock[i] + offsetOutside[i] ) inside = false;
                }
                
                if( inside ) {
                    
                    PointType GlobalCoordsBlock;
                    for( int i=0; i<GlobalCoordsBlock.size(); ++i){
                        GlobalCoordsBlock[i] = GlobalCoordsBlockOffset[i] - numDofsOffsetOutside * mesh.getMeshSize(i);
                    }
                    PointType GlobalCoordsSingleCell;
                    for( int i=0; i<GlobalCoordsSingleCell.size(); ++i){
                      int multDirection = static_cast<int> ( GlobalCoordsBlock[i] / lengthVec[i] );
                      GlobalCoordsSingleCell[i] = GlobalCoordsBlock[i] - multDirection * lengthVec[i];
                    }
                    int elementNumberSingleCell; PointType LocalCoordSingleCell;
                    conf.getLocalCoords ( GlobalCoordsSingleCell, elementNumberSingleCell, LocalCoordSingleCell );
                    materialBlockOutsideBone[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                    materialBlockOutsidePolymer[nodeIdxBlock] = discreteFctSingleCell.evaluate( mesh.getElement(elementNumberSingleCell), LocalCoordSingleCell );
                }
            }
            
            MeshWithData<MeshType> meshSaverBlockOutsideBone ( meshBlockOutside );
            meshSaverBlockOutsideBone.addScalarData ( materialBlockOutsideBone, "SolMaterialBlockOutsideBone", VERTEX_DATA );
            meshSaverBlockOutsideBone.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_BlockOutsideBone.vtk", saveDirectory.c_str() ), offsetOutside );
            
            MeshWithData<MeshType> meshSaverBlockOutsidePolymer ( meshBlockOutside );
            meshSaverBlockOutsidePolymer.addScalarData ( materialBlockOutsidePolymer, "SolMaterialBlockOutsidePolymer", VERTEX_DATA );
            meshSaverBlockOutsidePolymer.saveAsVTKSTRUCTUREDPOINTS ( aol::strprintf( "%s/SolMaterial_BlockOutsidePolymer.vtk", saveDirectory.c_str() ), offsetOutside );
        }
        
        //=====================================================================================================================
        //plot interface
        //=====================================================================================================================
        if(plotInterface){
            aol::consoleOutput( "plot interface" );
            parser.createSubDirectory("Interface");
            parser.createSubDirectory("Stress");
            
            //extract surface
            vtkInterfaceToSurface surfaceExtractor;
            const RealType threshold = parser.template get<RealType> ( "saving.thresholdInterface" );

            surfaceExtractor.getSurface( aol::strprintf( "%s/SolMaterial_OutsideBone.vtk", saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer.vtk", saveDirectory.c_str() ), threshold );
            surfaceExtractor.getSurface( aol::strprintf( "%s/SolMaterial_OutsidePolymer.vtk", saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", saveDirectory.c_str() ), threshold );
            
            surfaceExtractor.getSurface( aol::strprintf( "%s/SolMaterial_BlockOutsideBone.vtk", saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer.vtk", saveDirectory.c_str() ), threshold );
            surfaceExtractor.getSurface( aol::strprintf( "%s/SolMaterial_BlockOutsidePolymer.vtk", saveDirectory.c_str() ),
                                         aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", saveDirectory.c_str() ), threshold );
            
            //use loop subdivision for interfaces
            const int numSubdivLevels = 2;
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer.vtk", saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_InterfacePolymer_Subdiv%d.vtk", saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone_Subdiv%d.vtk", saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer.vtk", saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_BlockInterfacePolymer_Subdiv%d.vtk", saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
            surfaceExtractor.loopSubdivision( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", saveDirectory.c_str() ),
                                              aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone_Subdiv%d.vtk", saveDirectory.c_str(), numSubdivLevels ),
                                              numSubdivLevels );
                                              
            
            //! ========================================================================================
            //! plot von mises stresses on interface
            //! ========================================================================================
           
            aol::consoleOutput(  "plot von mises stresses" );
            MatOptConfigurator matOpConf ( parser, conf );
            const int numLoads = parser.template get<int> ( "AffineDisp.numLoads" );
            
            for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
                
                BonesInterface<MatOptConfigurator,BONE> interface ( matOpConf );
                
                VectorType displacementPeriodic ( conf.dimDomain * conf.getNumGlobalDofs() );
                cout << "size of displacementPeriodic = " << displacementPeriodic.size() << endl;
                aol::loadVectorFromFile<VectorType>( displacementPeriodic, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementBonePeriodic_Dir%d.txt", saveDirectory.c_str(),  loadIdx ) ); 
                
                VectorType displacementAffine ( conf.numAffineSymGradDofs );
                aol::loadVectorFromFile<VectorType>( displacementAffine, aol::strprintf( "%s/Deformation/SolMaterial_DisplacementBoneAffine_Dir%d.txt", saveDirectory.c_str(),  loadIdx ) ); 
                
                cout << "compute stress on interface" << endl;
                std::vector<RealType> stressVec, stressVecBlock;
                interface.evaluateStress( aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", saveDirectory.c_str() ).c_str(), material, displacementPeriodic, displacementAffine, stressVec, false  );
                interface.evaluateStress( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", saveDirectory.c_str() ).c_str(), material, displacementPeriodic, displacementAffine, stressVecBlock, true  );
                
                cout << "plot stress on interface" << endl;
                interface.plotStress( aol::strprintf( "%s/Interface/SolMaterial_InterfaceBone.vtk", saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Interface/SolMaterial_VonMises_InterfaceBone_Direction%d.vtk", saveDirectory.c_str(),  loadIdx ).c_str(),
                                      "VonMisesStresses", stressVec );
                
                interface.plotStress( aol::strprintf( "%s/Interface/SolMaterial_BlockInterfaceBone.vtk", saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Interface/SolMaterial_VonMises_BlockInterfaceBone_Direction%d.vtk", saveDirectory.c_str(),  loadIdx ).c_str(),
                                      "VonMisesStresses", stressVecBlock );
                
                //! ========================================================================================
                //plot stress on full 3d mesh
                //! ========================================================================================
                cout << "compute stress on 3d mesh" << endl;
                std::vector<RealType> stressVec3d, stressVecBlock3d;
                interface.evaluateStress( aol::strprintf( "%s/SolMaterial_Undeformed.vtk", saveDirectory.c_str() ).c_str(), material, displacementPeriodic, displacementAffine, stressVec3d, false  );
                interface.evaluateStress( aol::strprintf( "%s/SolMaterial_Block.vtk", saveDirectory.c_str() ).c_str(), material, displacementPeriodic, displacementAffine, stressVecBlock3d, true );
                cout << "plot stress on 3d mesh" << endl;
                interface.plotStress( aol::strprintf( "%s/SolMaterial_Undeformed.vtk", saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Stress/SolMaterial_VonMises_Bone_Direction%d.vtk", saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVec3d );
                
                interface.plotStress( aol::strprintf( "%s/SolMaterial_Block.vtk", saveDirectory.c_str() ).c_str(), 
                                      aol::strprintf( "%s/Stress/SolMaterial_VonMises_BlockBone_Direction%d.vtk", saveDirectory.c_str(), loadIdx ).c_str(),
                                      "VonMisesStresses", stressVecBlock3d );
                
            }
            
        }
        
        
        //=====================================================================================================================
        //compute homogenized tensor
        //=====================================================================================================================
        if(computeHomogenizedTensor){
            
            aol::consoleOutput( "compute homogenized tensor" );
            
            parser.createSubDirectory("HomogenizedTensor");
            
            const int dimDomain = ConfiguratorType::dimDomain;
            const int numAffineSymGradDofs = conf.numAffineSymGradDofs;
                       
            //define loads
            const int numLoadsAll = 6;
            std::vector<VectorType> affineDispBone, affineDispPolymer;
            for( int i=0; i < numLoadsAll; ++i ){
                VectorType affineDisp ( numAffineSymGradDofs ); affineDisp.setZero();
                affineDisp[i] = 0.1;                
                affineDispBone.push_back ( affineDisp ); affineDispPolymer.push_back ( affineDisp );
            }
            
            //!determine weightFunction for Multiple Load
            std::vector<RealType> weightFct_Bone_weightVec ( numLoadsAll ), weightFct_Polymer_weightVec ( numLoadsAll );
            for( int i=0; i<numLoadsAll; ++i){
              weightFct_Bone_weightVec[i] = 1.; weightFct_Polymer_weightVec[i] = 1.;   
            }
            
            RealType eps_area = 1.0, eps_factor = 1.0;
            MatOptConfigurator matOpConf ( parser, conf, eps_area, eps_factor, weightFct_Bone_weightVec, weightFct_Polymer_weightVec );
            
            //! solve
            OptimalDeformSolverMultipleLoad<MatOptConfigurator> OptDeformFinderAll( parser, matOpConf, material, quocHandler, affineDispBone, affineDispPolymer );

            ComplianceEnergyMultipleLoad_RegularizedMaxFunction<MatOptConfigurator> complianceOp ( OptDeformFinderAll );
            RealType complianceEnergy;
            complianceOp.apply( material, complianceEnergy );
            
            std::vector<RealType> homogenizedTensor_Bone ( numLoadsAll ), homogenizedTensor_Polymer ( numLoadsAll );
            for( int i=0; i<numLoadsAll; ++i ){
                homogenizedTensor_Bone[i] = complianceOp.getLastComplianceBone(i);
                homogenizedTensor_Polymer[i] = complianceOp.getLastCompliancePolymer(i);
            }
            
            
            //compute volume 
            shapeOptBonePolymerPeriodicHomogenization::VolumeConstraintPeriodicBC<MatOptConfigurator,BONE> volumeOpBone ( matOpConf, quocHandler );
            shapeOptBonePolymerPeriodicHomogenization::VolumeConstraintPeriodicBC<MatOptConfigurator,POLYMER> volumeOpPolymer ( matOpConf, quocHandler );
            RealType volumeBone; volumeOpBone.evaluateEnergy( material, volumeBone );
            RealType volumePolymer; volumeOpPolymer.evaluateEnergy( material, volumePolymer );
            
            MaterialOptimizationMultipleLoadEnergyInfo<MatOptConfigurator> energyInfoAll;
            energyInfoAll.setNumLoads( numLoadsAll );
            energyInfoAll.setComplianceBone( homogenizedTensor_Bone );
            energyInfoAll.setCompliancePolymer( homogenizedTensor_Polymer );
            energyInfoAll.setVolumeBone( volumeBone );
            energyInfoAll.setVolumePolymer( volumePolymer );
            
            TikzPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> tikzPlotterAll ( parser );
            tikzPlotterAll.template generateTexCodeForPaper<MatOptConfigurator>( "HomogenizedTensor/codeForPaper", energyInfoAll );

        }
        
        
        //=====================================================================================================================
        //Generate PDF with all results
        //=====================================================================================================================
        std::ofstream BashFilePDF ( aol::strprintf ( "%s/%s.sh", saveDirectory.c_str (), "generatePDF"  ) );
        BashFilePDF << "cd " << saveDirectory.c_str() << endl;
        BashFilePDF << "pdflatex Results" << ".tex" << endl;
        string systemCommand = "bash " + saveDirectory + "/generatePDF.sh";
        cout << "systemCommand = " << systemCommand << endl;
        bool failed;
        failed= ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
        if ( failed ) cerr << "programm returned an error." << endl;
    
    
    
  return 0;
}
