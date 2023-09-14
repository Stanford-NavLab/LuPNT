# Install headers
file(GLOB dynamics_headers "*.h")
install(FILES ${dynamics_headers} DESTINATION include/lupnt/dynamics)