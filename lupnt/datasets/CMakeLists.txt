# Install headers
file(GLOB numerics_headers "*.h")
install(FILES ${numerics_headers} DESTINATION include/lupnt/datasets)