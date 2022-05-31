# Add executables
add_executable(materialOptQuocLinElast quocMaterialOptimizationLinElast.cpp)

# Link the utility library
set(ELASTBODIESMATERIALOPTLINELASTLIBS ${SUITESPARSE_LIBRARIES} ${Boost_LIBRARIES})

target_link_libraries(materialOptQuocLinElast LINK_PUBLIC ${ELASTBODIESMATERIALOPTLINELASTLIBS})

if( USE_IPOPT )
 target_link_libraries(materialOptQuocLinElast LINK_PUBLIC ${IPOPT_LIBRARIES})
endif( USE_IPOPT )

if( USE_OPENMP )
 target_compile_definitions(materialOptQuocLinElast PUBLIC USE_OPENMP)
 set_target_properties(materialOptQuocLinElast PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_C_FLAGS})
endif( USE_OPENMP )

if( USE_VTK )
 target_link_libraries(materialOptQuocLinElast LINK_PUBLIC ${VTK_LIBRARIES} )
endif( USE_VTK )