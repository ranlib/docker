VER=v0.9.11

build:
	docker build -t dbest/pear:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/pear:$(VER) . 

push:
	docker push dbest/pear:$(VER)

clean:
	rm *~
