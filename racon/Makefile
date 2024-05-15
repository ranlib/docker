VER=v1.5.0
raconrounds:
	docker build -t dbest/raconrounds:latest -f racon.dockerfile .

raconrounds_publish:
	docker tag dbest/raconrounds:latest dbest/raconrounds:v2.2
	docker push dbest/raconrounds:v2.2

racon:
	docker build -t dbest/racon:$(VER) .

racon_no_cache:
	docker build --no-cache -t dbest/racon:$(VER) .

racon_push:
	docker push dbest/racon:$(VER)

clean:
	rm *~
