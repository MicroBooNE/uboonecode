set( ${product}_gdml_dir "${product}/${version}/gdml" )

# install gdml files
file(GLOB gdml_file *.gdml *.C *.xml README )
file(GLOB gdml_bin *.pl genmake )
install( FILES ${gdml_file} DESTINATION ${${product}_gdml_dir} )
install( PROGRAMS ${gdml_bin} DESTINATION ${${product}_gdml_dir} )


add_subdirectory(microboone)
add_subdirectory(GDMLSchema)
