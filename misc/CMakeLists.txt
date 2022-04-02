if (RSYNC)
  set (USER karney)
  set (DOCTOP $ENV{HOME}/web/geographiclib-web/htdocs)
  set (DOCROOT ${DOCTOP}/misc)
  set (WEBDEPLOY ${USER},geographiclib@web.sourceforge.net:./htdocs)

  add_custom_target (stage-misc
    COMMAND
      ${RSYNC} -av robots.txt google*.html index.html default.css ${DOCTOP}/
    COMMAND ${RSYNC} -av
      --exclude '.git*' --exclude '*~' --exclude '*.txt'
      --exclude 'google*.html' --exclude 'index.html' --exclude '*.png'
      --delete --delete-excluded . ${DOCROOT}/
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_custom_target (deploy-misc
    COMMAND
      ${RSYNC} --delete -av --exclude '*~' -e ssh ${DOCROOT} ${WEBDEPLOY}/ )

  add_custom_target (deploy-html
    COMMAND
      ${RSYNC} --delete -av --exclude '*~' -e ssh ${DOCTOP}/ ${WEBDEPLOY}/ )
endif ()
