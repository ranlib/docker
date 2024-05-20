VER=v1.8.3
build:
	docker build -t dbest/raven:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/raven:$(VER) . 

push:
	docker push dbest/raven:$(VER)

clean:
	rm *~
