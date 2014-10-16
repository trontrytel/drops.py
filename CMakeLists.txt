cmake_minimum_required(VERSION 2.8)
project(drops.py NONE)

find_package(PythonInterp REQUIRED)

install(
  FILES drops.py 
  DESTINATION bin 
  PERMISSIONS 
    OWNER_EXECUTE OWNER_WRITE OWNER_READ 
    GROUP_EXECUTE GROUP_READ
    WORLD_EXECUTE WORLD_READ
)

execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
    "import site, sys; sys.stdout.write(site.PREFIXES[-1])"
   OUTPUT_VARIABLE PYTHON_PREFIX
)
file ( TO_CMAKE_PATH "${PYTHON_PREFIX}" PYTHON_PREFIX )
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
    "import site, sys; sys.stdout.write(site.getsitepackages()[-1])"
   OUTPUT_VARIABLE PYTHON_SITE_DIR
)
file ( TO_CMAKE_PATH "${PYTHON_SITE_DIR}" PYTHON_SITE_DIR )
string ( REGEX REPLACE "^${PYTHON_PREFIX}/" ""
   PYTHON_SITE_DIR "${PYTHON_SITE_DIR}"
)

install( 
  DIRECTORY drops_py/
  DESTINATION "${PYTHON_SITE_DIR}/drops_py"
  PATTERN "drops_py/*.py"
)
