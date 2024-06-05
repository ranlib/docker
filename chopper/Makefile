VER=v0.8.0

build:
	docker build -t dbest/chopper:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/chopper:$(VER) . 

push:
	docker push dbest/chopper:$(VER)
