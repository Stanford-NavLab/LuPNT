# Install headers
file(GLOB physics_headers "*.h")
install(FILES ${physics_headers} DESTINATION include/lupnt/physics)