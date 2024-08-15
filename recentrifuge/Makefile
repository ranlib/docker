VER=v1.14.1

build:
	docker build -t dbest/recentrifuge:$(VER) -f Dockerfile .

build_no_cache:
	docker build --no-cache -t dbest/recentrifuge:$(VER) -f Dockerfile .

push:
	docker push  dbest/recentrifuge:$(VER)

clean:
	rm *~
