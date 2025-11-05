VERSION=v3.0.0
TAG = dbest/python3:$(VERSION)

build:
	docker build -t $(TAG) .

build_no_cache:
	docker build --no-cache -t $(TAG) .

push:
	docker push $(TAG)
