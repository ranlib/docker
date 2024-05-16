IMAGE_NAME = asm
VERSION = v0.1.13

TAG = dbest/$(IMAGE_NAME):$(VERSION)

all: | build push

build:
	docker build -t $(TAG) .

build_no_cache:
	docker build --no-cache -t $(TAG) .

push:
	docker push $(TAG)
