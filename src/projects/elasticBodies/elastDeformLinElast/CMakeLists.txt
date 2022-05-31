# Add executables
add_executable(optimalDeformQuocLinElast quocOptimalDeformLinElast.cpp)

# Link the utility library
set(ELASTBODIESDEFORMLINELASTLIBS ${SUITESPARSE_LIBRARIES} ${Boost_LIBRARIES})

target_link_libraries(optimalDeformQuocLinElast LINK_PUBLIC ${ELASTBODIESDEFORMLINELASTLIBS})
if( USE_OPENMP )
 target_compile_definitions(optimalDeformQuocLinElast PUBLIC USE_OPENMP)
 set_target_properties(optimalDeformQuocLinElast PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
endif( USE_OPENMP )

if( USE_VTK )
 target_link_libraries(optimalDeformQuocLinElast LINK_PUBLIC ${VTK_LIBRARIES} )
endif( USE_VTK )