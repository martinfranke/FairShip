# Create a library called "libPassive" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/passive
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}

)
 
link_directories( ${LINK_DIRECTORIES})

set(SRCS
ShipCave.cxx            
ShipGeoCave.cxx        
ShipMagnet.cxx         
ShipChamber.cxx
ShipTargetStation.cxx
ShipMuonShield.cxx
ShipPassiveContFact.cxx
ShipTAUMagneticSpectrometer.cxx
ShipGoliath.cxx
)

Set(HEADERS )
Set(LINKDEF PassiveLinkDef.h)
Set(LIBRARY_NAME Passive)
Set(DEPENDENCIES Base GeoBase ParBase Geom Cint Core)

GENERATE_LIBRARY()
