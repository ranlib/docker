VER=v2.2.1

build:
	docker build -t dbest/ncbi-scrubber:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/ncbi-scrubber:$(VER) .

push:
	docker push dbest/ncbi-scrubber:$(VER)
