set (SOURCES
  oindex.html
  LICENSE.txt
  robots.txt
  google16ac598d95b3d4df.html
  default.css
  geod.html
  tm.html
  geod-addenda.html
  bessel-errata.html
  jacobi-errata.html
  tm-addenda.html
  tm-figs.pdf
  tm-grid.kmz)

if (RSYNC)
  set (USER karney)
  set (DOCROOT $ENV{HOME}/web/geographiclib-web/htdocs)
  set (WEBDEPLOY ${USER},geographiclib@web.sourceforge.net:./htdocs)

  add_custom_target (stage-misc
    COMMAND ${RSYNC} -av ${SOURCES} ${DOCROOT}/
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_custom_target (deploy-html
    COMMAND
      ${RSYNC} --delete -av --exclude '*~' -e ssh ${DOCROOT}/ ${WEBDEPLOY}/ )
endif ()
