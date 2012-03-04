set(grappa_dist_src 
	invdist.c 
	med_util.c
	uf.c
)

set(grappa_capr_src 
	inversion_median_alberto.c
)

set_source_files_properties(inversion_median_alberto.c
	PROPERTIES
	LANGUAGE CXX )


add_library (grappa_dist SHARED ${grappa_dist_src})
add_library (grappa_capr SHARED ${grappa_capr_src})

#add_library (caprara SHARED caprara.cpp)

set_target_properties(grappa_capr 
	PROPERTIES 
	LINKER_LANGUAGE CXX 
	LANGUAGE CXX)

ADD_EXECUTABLE (caprara_main caprara.cpp)

# link the executable to the library
target_link_libraries (caprara_main grappa_capr grappa_dist)

