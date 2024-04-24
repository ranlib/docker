VER = v0.3.0

build:
	docker build -t dbest/dorado:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/dorado:$(VER) .

push:
	docker push dbest/dorado:$(VER)
