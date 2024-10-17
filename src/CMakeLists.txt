# BUILD LIBRARY 
file(GLOB SRC "*.cpp" ) 
include_directories($ENV{ITERATEDOKT}/src)

ROOT_GENERATE_DICTIONARY(G__ITERATEDOKT
    LINKDEF ITERATEDOKTLinkDef.h
)
add_library(ITERATEDOKT SHARED ${SRC} G__ITERATEDOKT.cxx )
target_link_libraries( ITERATEDOKT ROOT::Core ROOT::MathMore )
target_link_libraries( ITERATEDOKT Boost::boost)

install( TARGETS ITERATEDOKT
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libITERATEDOKT_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

# BUILD EXECUTABLE
add_executable( iteratedOKT "cling/iteratedOKT.cpp" )
target_link_libraries( iteratedOKT ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( iteratedOKT ITERATEDOKT )