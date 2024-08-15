VER=v1.0.4.1

build:
	docker build -t dbest/centrifuge:$(VER) .

build_no_cache:
	docker build --no-cache -t dbest/centrifuge:$(VER) .

push:
	docker push dbest/centrifuge:$(VER)

clean:
	rm *~

