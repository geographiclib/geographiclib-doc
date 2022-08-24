add_custom_command (OUTPUT biblio.html
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/biblio.sh
    ${CMAKE_CURRENT_SOURCE_DIR}/geodesic-biblio.txt
  COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/biblio.sh
    < ${CMAKE_CURRENT_SOURCE_DIR}/geodesic-biblio.txt > biblio.html )

add_custom_target (bib DEPENDS biblio.html)

if (RSYNC)
  set (USER karney)
  set (DOCROOT $ENV{HOME}/web/geographiclib-web/htdocs/geodesic-papers)
  set (WEBDEPLOY ${USER},geographiclib@web.sourceforge.net:./htdocs)

  add_custom_target (stage-bib
    COMMAND ${RSYNC} -av --delete --exclude biblio.html --exclude '*.txt'
      --exclude biblio.sh --exclude '.git*' --exclude '*~' . ${DOCROOT}/
    COMMAND ${RSYNC} -av ${CMAKE_CURRENT_BINARY_DIR}/biblio.html ${DOCROOT}/
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_dependencies (stage-bib bib)

  add_custom_target (deploy-bib
    COMMAND
      ${RSYNC} --delete -av --exclude '*~' -e ssh ${DOCROOT} ${WEBDEPLOY}/ )
endif ()
