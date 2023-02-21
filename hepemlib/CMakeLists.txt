add_library(hepemlib
    include/hepemlib/PhysicsList.hh
    include/hepemlib/PhysListHepEm.hh
    include/hepemlib/PhysListHepEmTracking.hh
    src/PhysicsList.cc
    src/PhysListHepEm.cc
    src/PhysListHepEmTracking.cc
)

target_include_directories(hepemlib PUBLIC include)
target_link_libraries(hepemlib ${Geant4_LIBRARIES} ${G4HepEm_LIBRARIES})
