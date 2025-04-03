VER = v2.5.2
build:
	docker build -t dbest/nextdenovo:$(VER) .

build_nocache:
	docker build --no-cache -t dbest/nextdenovo:$(VER) .

push:
	docker push dbest/nextdenovo:$(VER)
