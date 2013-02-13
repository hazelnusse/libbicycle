set(MATHJAX_RELPATH "http://cdn.mathjax.org/mathjax/latest")
if (NOT BICYCLE_DOXYGEN_MATHJAX_USECDN)
  if (NOT BICYCLE_DOXYGEN_MATHJAX_LOCAL STREQUAL "")
    set(MATHJAX_RELPATH ${BICYCLE_DOXYGEN_MATHJAX_LOCAL})
    message(${MATHJAX_RELPATH})
  endif()
endif()
# prepare doxygen configuration file
configure_file(${BICYCLE_SOURCE_DIR}/doc/Doxyfile.in
               ${BICYCLE_BINARY_DIR}/doc/Doxyfile)

configure_file(${BICYCLE_SOURCE_DIR}/doc/DoxygenLayout.xml.in
               ${BICYCLE_BINARY_DIR}/doc/DoxygenLayout.xml)

configure_file(${BICYCLE_SOURCE_DIR}/doc/bicycle_header.html.in
               ${BICYCLE_BINARY_DIR}/doc/bicycle_header.html)

configure_file(${BICYCLE_SOURCE_DIR}/doc/bicycle_footer.html.in
               ${BICYCLE_BINARY_DIR}/doc/bicycle_footer.html)

configure_file(${BICYCLE_SOURCE_DIR}/doc/bicycle_stylesheet.css.in
               ${BICYCLE_BINARY_DIR}/doc/bicycle_stylesheet.css)

add_custom_target(doc
                  ${DOXYGEN_EXECUTABLE}
                  ${BICYCLE_BINARY_DIR}/doc/Doxyfile
                  WORKING_DIRECTORY ${BICYCLE_BINARY_DIR}
                  COMMENT "Generating API Documentation with Doxygen"
                  VERBATIM)
