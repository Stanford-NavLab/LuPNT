# Install headers
file(GLOB agents_headers "*.h")
install(FILES ${agents_headers} DESTINATION include/lupnt/agents)