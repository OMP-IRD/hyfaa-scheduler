docker run --rm -it \
  -e HYDROWEB_USER="remi.jugier@magellium.fr" \
  -e HYDROWEB_PASSWORD="hydroMGB31"  \
  -v $(pwd)/work_configurations/operational_niger_gsmap:/work \
  hyfaa:4.0 $@
