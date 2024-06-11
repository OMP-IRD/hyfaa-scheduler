docker run --rm -it \
  -e EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY="yourverysecretapikey" \
  -v $(pwd)/work_configurations/operational_niger_gsmap:/work \
  hyfaa:4.0 $@
