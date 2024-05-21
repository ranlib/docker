VER=v0.3.1
build:
	docker build -t dbest/polca:$(VER) . 

build_no_cache:
	docker build --no-cache -t dbest/polca:$(VER) . 

push:
	docker push dbest/polca:$(VER)

clean:
	rm *~
