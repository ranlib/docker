VER=v4.1.1

build:
	docker build -t dbest/metaphlan:$(VER) -f Dockerfile .

build_no_cache:
	docker build --no-cache -t dbest/metaphlan:$(VER) -f Dockerfile .

push:
	docker push  dbest/metaphlan:$(VER)

clean:
	rm *~
