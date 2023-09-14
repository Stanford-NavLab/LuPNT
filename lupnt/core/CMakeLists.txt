# Install headers
file(GLOB core_headers "*.h")
install(FILES ${core_headers} DESTINATION include/lupnt/core)