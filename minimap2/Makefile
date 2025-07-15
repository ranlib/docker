VERSION = 2.30
TAG = dbest/minimap2:$(VERSION)

build:
	docker build -t $(TAG) .

build_no_cache:
	docker build --no-cache -t $(TAG) .

push:
	docker push $(TAG)
