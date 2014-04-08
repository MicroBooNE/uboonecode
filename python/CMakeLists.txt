FILE(GLOB python_modules *.py)
install(FILES ${python_modules} DESTINATION ${product}/${version}/python)
