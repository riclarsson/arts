set(files
  guide.rst
  guide.theory.rst
  guide.theory.absorption.rst
  guide.theory.absorption.line-by-line.rst
  guide.theory.absorption.spectral_radiance.rst
  guide.theory.absorption.propagation_matrix.rst
  guide.user.rst
  develop.rst
)

foreach(fil ${files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${fil} ${CMAKE_BINARY_DIR}/python/doc/source/${fil} COPYONLY)
endforeach(fil ${files})
