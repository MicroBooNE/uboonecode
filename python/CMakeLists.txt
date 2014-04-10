# All of the .py files in this directory are pyton modules.
# We distinguish executable and non-executable .py files.
# The former are isntalled with execute permission.
# The latter are installed without execute permission.

# Executable python files.

LIST(APPEND exes extractor_dict.py)

# Non-executable python files.

LIST(APPEND nonexes uboone_utilities.py) 

message(STATUS "Executable python modules ${exes}")
message(STATUS "Non-executable python modules ${nonexes}")

install(PROGRAMS ${exes} DESTINATION ${product}/${version}/python)
install(FILES ${nonexes} DESTINATION ${product}/${version}/python)
install_scripts(LIST ${exes})
