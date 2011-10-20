#PROJECT(direct Fortran)

enable_language(Fortran)

ADD_LIBRARY(DIRect 
  Direct.f 
  DIRserial.f 
  DIRsubrout.f 
  DirectHead.f 
)

INSTALL(TARGETS DIRect
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)