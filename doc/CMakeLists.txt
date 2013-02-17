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

configure_file(${BICYCLE_SOURCE_DIR}/doc/html/DoxygenLayout.xml.in
               ${BICYCLE_BINARY_DIR}/doc/html/DoxygenLayout.xml)

configure_file(${BICYCLE_SOURCE_DIR}/doc/html/header.html.in
               ${BICYCLE_BINARY_DIR}/doc/html/header.html)

configure_file(${BICYCLE_SOURCE_DIR}/doc/html/footer.html.in
               ${BICYCLE_BINARY_DIR}/doc/html/footer.html)

configure_file(${BICYCLE_SOURCE_DIR}/doc/html/stylesheet.css.in
               ${BICYCLE_BINARY_DIR}/doc/html/stylesheet.css)

configure_file(${BICYCLE_SOURCE_DIR}/doc/latex/header.tex.in
               ${BICYCLE_BINARY_DIR}/doc/latex/header.tex)

configure_file(${BICYCLE_SOURCE_DIR}/doc/latex/footer.tex.in
               ${BICYCLE_BINARY_DIR}/doc/latex/footer.tex)

configure_file(${BICYCLE_SOURCE_DIR}/doc/latex/doxygen.sty.in
               ${BICYCLE_BINARY_DIR}/doc/latex/doxygen.sty)


add_custom_target(doc
                  ${DOXYGEN_EXECUTABLE}
                  ${BICYCLE_BINARY_DIR}/doc/Doxyfile
                  WORKING_DIRECTORY ${BICYCLE_BINARY_DIR}
                  COMMENT "Generating API Documentation with Doxygen"
                  VERBATIM)
