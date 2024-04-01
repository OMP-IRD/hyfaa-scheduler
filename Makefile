SHELL := /bin/bash
TAG=latest
IMAGE=pigeosolutions/hyfaa-scheduler
REV=`git rev-parse --short HEAD`
DATE=`date +%Y%m%d-%H%M`

all: docker-build docker-push

docker-build:
	docker build -f Dockerfile-pipinstall -t ${IMAGE}:latest . ;\
	echo built ${IMAGE}:latest ;\
	docker tag  ${IMAGE}:latest ${IMAGE}:${DATE}-${REV} ;\
	echo built ${IMAGE}:${DATE}-${REV}

docker-push:
	docker push ${IMAGE}:latest ;\
	echo pushed ${IMAGE}:latest ;\
	docker push ${IMAGE}:${DATE}-${REV} ;\
	echo pushed ${IMAGE}:${DATE}-${REV}
