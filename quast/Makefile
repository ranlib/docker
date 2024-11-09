VER = v5.2.0
build:
	docker build -t dbest/quast:$(VER) .

build_nocache:
	docker build --no-cache -t dbest/quast:$(VER) .

push:
	docker push dbest/quast:$(VER)
