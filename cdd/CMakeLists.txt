include_directories(${GMP_INCLUDE_PATH})

set(SOURCE_FILES

#### Header ####

cdd.h
cddmp.h
cddtypes.h

cdd_f.h
cddmp_f.h
cddtypes_f.h

setoper.h

#### Sources ####

cddcore.c
cddio.c
cddlib.c
cddlp.c
cddmp.c
cddproj.c

cddcore_f.c
cddio_f.c
cddlib_f.c
cddlp_f.c
cddmp_f.c
cddproj_f.c

setoper.c

)

add_definitions(-DGMPRATIONAL)

add_library(cdd STATIC ${SOURCE_FILES})