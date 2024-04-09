set (SOURCEDIR "${CMAKE_CURRENT_SOURCE_DIR}")
set (DOCSOURCES
  ${SOURCEDIR}/index.rst
  ${SOURCEDIR}/doc/library.rst
  ${SOURCEDIR}/doc/research.rst
  ${SOURCEDIR}/doc/triaxial.rst)

add_custom_target (doc ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
  DEPENDS ${DOCSOURCES} ${SOURCEDIR}/conf.py
  COMMAND ${CMAKE_COMMAND} -E rm -rf html
  COMMAND ${SPHINX} -v -b html ${SOURCEDIR} html > docs.log
  COMMENT "Generating html documentation tree")

if (RSYNC)
  set (USER karney)
  set (DOCROOT $ENV{HOME}/web/geographiclib-web/htdocs)
  set (WEBDEPLOY ${USER},geographiclib@web.sourceforge.net:./htdocs)
  set (SPHINX_OUTPUT
    doc
    genindex.html
    index.html
    objects.inv
    search.html
    searchindex.js
    _sources
    _static
    )

  add_custom_target (stage-doc
    COMMAND ${RSYNC} -av --delete ${SPHINX_OUTPUT} ${DOCROOT}/
    WORKING_DIRECTORY html)
  add_dependencies (stage-doc doc)

  add_custom_target (deploy-doc
    COMMAND
      ${RSYNC} --delete -av --exclude '*~' -e ssh ${SPHINX_OUTPUT} ${WEBDEPLOY}/
    WORKING_DIRECTORY ${DOCROOT})

endif ()
