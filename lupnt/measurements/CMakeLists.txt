# Install headers
file(GLOB gnss_headers "*.h")
install(FILES ${gnss_headers} DESTINATION include/lupnt/gnss)